#include <Rcpp.h>

using namespace Rcpp;

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "gauss.h"
#include "snp.h"
#include "util.h"

namespace {

struct SuperPopLDWindowCleanup {
  std::map<MapKey, Snp*, LessThanMapKey>& snp_map;
  std::vector<Snp*>& snp_vec;
  bool genotype_loaded;

  SuperPopLDWindowCleanup(std::map<MapKey, Snp*, LessThanMapKey>& map_ref,
                          std::vector<Snp*>& vec_ref)
      : snp_map(map_ref), snp_vec(vec_ref), genotype_loaded(false) {}

  ~SuperPopLDWindowCleanup() {
    if (genotype_loaded) {
      FreeGenotype(snp_vec);
    }
    for (std::map<MapKey, Snp*, LessThanMapKey>::iterator it = snp_map.begin(); it != snp_map.end();) {
      (it->second)->ClearSnp();
      delete it->second;
      snp_map.erase(it++);
    }
  }
};

// Resolves an optional `sample_keep` list (names = population abbreviations,
// values = 1-based positions of individuals within that population's genotype
// string) into 0-based positions per selected population, in the selected-pop
// order ReadGenotype uses. Returns an empty vector when `sample_keep` is NULL.
// Every selected population must be named, so a typo can't silently fall back
// to using all of a population's samples.
std::vector<std::vector<int> > resolve_sample_keep(const Arguments& args,
                                                   const Rcpp::Nullable<Rcpp::List>& sample_keep) {
  std::vector<std::vector<int> > keep;
  if (sample_keep.isNull()) {
    return keep;
  }
  Rcpp::List sk(sample_keep);
  if (sk.size() > 0 && Rf_isNull(sk.names())) {
    Rcpp::stop("'sample_keep' must be a named list keyed by population abbreviation.");
  }
  Rcpp::CharacterVector sk_names = sk.size() > 0 ? Rcpp::CharacterVector(sk.names()) : Rcpp::CharacterVector(0);

  for (int i = 0; i < args.num_pops; ++i) {
    if (!args.pop_flag_vec[i]) continue;
    const std::string& pop = args.ref_pop_vec[i];
    int found = -1;
    for (int j = 0; j < sk_names.size(); ++j) {
      if (Rcpp::as<std::string>(sk_names[j]) == pop) {
        if (found >= 0) Rcpp::stop("'sample_keep' names population " + pop + " more than once.");
        found = j;
      }
    }
    if (found < 0) {
      Rcpp::stop("'sample_keep' has no entry for selected population " + pop + ".");
    }
    Rcpp::IntegerVector pos(sk[found]);
    const int pop_n = args.ref_pop_size_vec[i];
    std::vector<int> pos0;
    pos0.reserve(pos.size());
    for (int k = 0; k < pos.size(); ++k) {
      if (pos[k] == NA_INTEGER || pos[k] < 1 || pos[k] > pop_n) {
        Rcpp::stop("'sample_keep' position out of range for population " + pop + " (size " +
                   std::to_string(pop_n) + ").");
      }
      pos0.push_back(pos[k] - 1);
    }
    std::vector<int> sorted = pos0;
    std::sort(sorted.begin(), sorted.end());
    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) {
      Rcpp::stop("'sample_keep' repeats a position for population " + pop + ".");
    }
    if (pos0.empty()) {
      Rcpp::stop("'sample_keep' keeps no individuals for population " + pop + ".");
    }
    keep.push_back(pos0);
  }
  return keep;
}

// Reduces each SNP's per-population genotype strings to the kept positions.
void subset_genotypes(std::vector<Snp*>& snp_vec, const std::vector<std::vector<int> >& keep) {
  for (size_t s = 0; s < snp_vec.size(); ++s) {
    std::vector<std::string>& geno_vec = snp_vec[s]->GetGenotypeVec();
    if (geno_vec.size() != keep.size()) {
      Rcpp::stop("Genotype population count does not match 'sample_keep' for SNP " + snp_vec[s]->GetRsid() + ".");
    }
    for (size_t p = 0; p < geno_vec.size(); ++p) {
      const std::string& full = geno_vec[p];
      std::string sub;
      sub.reserve(keep[p].size());
      for (size_t k = 0; k < keep[p].size(); ++k) {
        if (static_cast<size_t>(keep[p][k]) >= full.size()) {
          Rcpp::stop("Genotype string shorter than population size for SNP " + snp_vec[s]->GetRsid() + ".");
        }
        sub.push_back(full[keep[p][k]]);
      }
      geno_vec[p] = sub;
    }
  }
}

// Shared setup for one superpopulation window: selects the populations,
// reads the index (restricted to `keep_bp` if given), loads genotypes, applies
// `sample_keep`, and sets args.pop_wgt_vec to sample-size weights of the
// individuals actually kept. Returns the number of selected populations;
// `selected_pops` receives their names in genotype order.
int load_superpop_window(Arguments& args,
                         std::map<MapKey, Snp*, LessThanMapKey>& snp_map,
                         std::vector<Snp*>& snp_vec,
                         SuperPopLDWindowCleanup& cleanup,
                         const Rcpp::Nullable<Rcpp::IntegerVector>& keep_bp,
                         const Rcpp::Nullable<Rcpp::List>& sample_keep,
                         std::vector<std::string>& selected_pops) {
  read_ref_desc(args);

  // Existing GAUSS helper supports both population and superpopulation labels.
  init_pop_flag_vec(args);

  std::vector<double> pop_n;
  for (int i = 0; i < args.num_pops; ++i) {
    if (args.pop_flag_vec[i]) {
      selected_pops.push_back(args.ref_pop_vec[i]);
      pop_n.push_back(static_cast<double>(args.ref_pop_size_vec[i]));
    }
  }
  const int n_selected_pops = static_cast<int>(selected_pops.size());

  const std::vector<std::vector<int> > keep_samples = resolve_sample_keep(args, sample_keep);
  if (!keep_samples.empty()) {
    for (int p = 0; p < n_selected_pops; ++p) {
      pop_n[p] = static_cast<double>(keep_samples[p].size());
    }
  }

  double selected_sample_total = 0.0;
  for (int p = 0; p < n_selected_pops; ++p) selected_sample_total += pop_n[p];

  if (n_selected_pops == 0 || selected_sample_total <= 0.0) {
    Rcpp::stop("Invalid superpopulation label: no reference populations selected.");
  }

  // Population weights in the same selected-pop order used by ReadGenotype,
  // proportional to the individuals actually used.
  args.pop_wgt_vec.clear();
  for (int p = 0; p < n_selected_pops; ++p) {
    args.pop_wgt_vec.push_back(pop_n[p] / selected_sample_total);
  }

  ReadReferenceIndex(snp_map, args);

  // SNPs left out here stay in snp_map, so the cleanup still deletes them; they
  // never have genotypes loaded, so leaving them out of snp_vec is safe.
  std::vector<long long int> keep;
  const bool filter_bp = keep_bp.isNotNull();
  if (filter_bp) {
    Rcpp::IntegerVector kb(keep_bp);
    keep.assign(kb.begin(), kb.end());
    std::sort(keep.begin(), keep.end());
  }
  for (std::map<MapKey, Snp*, LessThanMapKey>::iterator it = snp_map.begin(); it != snp_map.end(); ++it) {
    if (filter_bp && !std::binary_search(keep.begin(), keep.end(), static_cast<long long int>(it->second->GetBp()))) {
      continue;
    }
    snp_vec.push_back(it->second);
  }

  if (snp_vec.empty()) {
    Rcpp::stop("No SNPs found in requested region for provided reference index.");
  }

  ReadGenotype(snp_vec, args);
  cleanup.genotype_loaded = true;

  if (!keep_samples.empty()) {
    subset_genotypes(snp_vec, keep_samples);
  }
  return n_selected_pops;
}

Arguments make_window_args(int chr, long long int start_bp, long long int end_bp,
                           std::string superpopulation_label,
                           const std::string& reference_index_file,
                           const std::string& reference_data_file,
                           const std::string& reference_pop_desc_file) {
  if (start_bp > end_bp) {
    Rcpp::stop("Invalid region: 'start_bp' must be <= 'end_bp'.");
  }
  std::transform(superpopulation_label.begin(), superpopulation_label.end(),
                 superpopulation_label.begin(),
                 [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

  Arguments args;
  args.chr = chr;
  args.start_bp = start_bp;
  args.end_bp = end_bp;
  args.wing_size = 0;
  args.study_pop = superpopulation_label;
  args.reference_index_file = reference_index_file;
  args.reference_data_file = reference_data_file;
  args.reference_pop_desc_file = reference_pop_desc_file;
  return args;
}

}  // namespace

//' Compute superpopulation-specific LD for one reference window
//'
//' Computes a signed LD correlation matrix for one genomic region and one selected
//' superpopulation (AFR/AMR/ASN/EUR/SAS) directly from GAUSS 33KG reference files,
//' without GWAS summary statistics input.
//'
//' @param chr Chromosome number.
//' @param start_bp Start base-pair position of the region.
//' @param end_bp End base-pair position of the region.
//' @param superpopulation_label Superpopulation label (e.g. AFR/AMR/ASN/EUR/SAS).
//' @param reference_index_file Reference panel index file.
//' @param reference_data_file Reference panel genotype data file.
//' @param reference_pop_desc_file Reference panel population description file.
//' @param maf_cutoff Optional MAF threshold. Defaults to 0.01.
//' @param missing_cutoff Optional missingness threshold in [0,1]. Defaults to 1.0.
//' @param keep_bp Optional integer vector of base-pair positions. When given,
//'   only reference SNPs at these positions are read and used, which avoids
//'   reading genotypes and computing LD for SNPs that will be discarded (e.g.
//'   when restricting a panel to HapMap3 SNPs). Positions only: allele
//'   matching at multi-allelic sites is left to the caller.
//' @param min_snps Minimum number of SNPs that must survive filtering, else an
//'   "insufficient SNPs" error is raised. Defaults to 2; 1 allows single-SNP
//'   (1x1) blocks.
//' @param sample_keep Optional named list restricting which individuals are
//'   used: names are the superpopulation's population abbreviations (every
//'   one must appear), values are 1-based positions of individuals within that
//'   population's genotype string. Population weights then follow the number
//'   of individuals kept. Used to build panels from one part of a sample split
//'   (see \code{\link{split_reference_samples}}).
//' @return A list with `snplist` (data frame: rsid, chr, bp, a1, a2, af1pop)
//'   and `cormat` (signed correlation matrix).
// [[Rcpp::export]]
Rcpp::List computeSuperPopLDWindow(
    int chr,
    long long int start_bp,
    long long int end_bp,
    std::string superpopulation_label,
    std::string reference_index_file,
    std::string reference_data_file,
    std::string reference_pop_desc_file,
    Rcpp::Nullable<double> maf_cutoff = R_NilValue,
    Rcpp::Nullable<double> missing_cutoff = R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> keep_bp = R_NilValue,
    int min_snps = 2,
    Rcpp::Nullable<Rcpp::List> sample_keep = R_NilValue) {

  if (min_snps < 1) {
    Rcpp::stop("Invalid 'min_snps': must be at least 1.");
  }
  const size_t min_kept = static_cast<size_t>(min_snps);

  const double maf_min = maf_cutoff.isNotNull() ? Rcpp::as<double>(maf_cutoff) : 0.01;
  const double missing_max = missing_cutoff.isNotNull() ? Rcpp::as<double>(missing_cutoff) : 1.0;

  if (!std::isfinite(maf_min) || maf_min < 0.0 || maf_min >= 0.5) {
    Rcpp::stop("Invalid 'maf_cutoff': must be finite and in [0, 0.5)." );
  }
  if (!std::isfinite(missing_max) || missing_max < 0.0 || missing_max > 1.0) {
    Rcpp::stop("Invalid 'missing_cutoff': must be finite and in [0, 1].");
  }

  Arguments args = make_window_args(chr, start_bp, end_bp, superpopulation_label,
                                    reference_index_file, reference_data_file,
                                    reference_pop_desc_file);

  std::map<MapKey, Snp*, LessThanMapKey> snp_map;
  std::vector<Snp*> snp_vec;
  SuperPopLDWindowCleanup cleanup(snp_map, snp_vec);
  std::vector<std::string> selected_pops;
  const int n_selected_pops = load_superpop_window(args, snp_map, snp_vec, cleanup,
                                                   keep_bp, sample_keep, selected_pops);

  std::vector<Snp*> qc_snps;
  std::vector<double> qc_af1;
  qc_snps.reserve(snp_vec.size());
  qc_af1.reserve(snp_vec.size());

  for (size_t i = 0; i < snp_vec.size(); ++i) {
    Snp* snp = snp_vec[i];
    std::vector<std::string>& geno_vec = snp->GetGenotypeVec();

    if (static_cast<int>(geno_vec.size()) != n_selected_pops) {
      Rcpp::stop("Superpopulation mode expected genotype_vec size " + std::to_string(n_selected_pops) +
                 ", but got " + std::to_string(geno_vec.size()) + " for SNP " + snp->GetRsid() + ".");
    }

    int total_n = 0;
    int nonmissing_n = 0;
    int missing_n = 0;
    double dosage_sum = 0.0;

    for (size_t p = 0; p < geno_vec.size(); ++p) {
      const std::string& geno = geno_vec[p];
      total_n += static_cast<int>(geno.size());
      for (size_t j = 0; j < geno.size(); ++j) {
        char c = geno[j];
        if (c == '0' || c == '1' || c == '2') {
          dosage_sum += static_cast<double>(c - '0');
          nonmissing_n++;
        } else {
          missing_n++;
        }
      }
    }

    if (total_n == 0 || nonmissing_n == 0) {
      continue;
    }

    const double missing_frac = static_cast<double>(missing_n) / static_cast<double>(total_n);
    const double af1pop = dosage_sum / (2.0 * static_cast<double>(nonmissing_n));
    const double maf = std::min(af1pop, 1.0 - af1pop);

    if (missing_frac <= missing_max && maf >= maf_min) {
      qc_snps.push_back(snp);
      qc_af1.push_back(af1pop);
    }
  }

  if (qc_snps.size() < min_kept) {
    Rcpp::stop("Fewer than " + std::to_string(min_snps) + " SNPs remain after missingness/MAF filtering.");
  }

  std::vector<Snp*> kept_snps;
  std::vector<double> kept_af1;
  std::vector<double> snp_sd;
  kept_snps.reserve(qc_snps.size());
  kept_af1.reserve(qc_snps.size());
  snp_sd.reserve(qc_snps.size());

  for (size_t i = 0; i < qc_snps.size(); ++i) {
    Snp* snp = qc_snps[i];
    std::vector<std::string>& geno_vec = snp->GetGenotypeVec();

    // CalWgtCov expects dosage digits only.
    bool has_missing = false;
    for (size_t p = 0; p < geno_vec.size() && !has_missing; ++p) {
      const std::string& geno = geno_vec[p];
      for (size_t j = 0; j < geno.size(); ++j) {
        char c = geno[j];
        if (!(c == '0' || c == '1' || c == '2')) {
          has_missing = true;
          break;
        }
      }
    }

    if (has_missing) {
      continue;
    }

    const double v = CalWgtCov(geno_vec, geno_vec, args.pop_wgt_vec);
    if (!std::isfinite(v) || v <= 0.0) {
      continue;
    }

    kept_snps.push_back(snp);
    kept_af1.push_back(qc_af1[i]);
    snp_sd.push_back(std::sqrt(v));
  }

  if (kept_snps.size() < min_kept) {
    Rcpp::stop("Fewer than " + std::to_string(min_snps) + " SNPs remain after variance/non-missing genotype checks.");
  }

  const int m = static_cast<int>(kept_snps.size());
  Rcpp::NumericMatrix cor_mat(m, m);

  for (int i = 0; i < m; ++i) {
    cor_mat(i, i) = 1.0;
    for (int j = i + 1; j < m; ++j) {
      const double cov = CalWgtCov(kept_snps[i]->GetGenotypeVec(), kept_snps[j]->GetGenotypeVec(), args.pop_wgt_vec);
      const double cor = cov / (snp_sd[i] * snp_sd[j]);
      if (!std::isfinite(cor)) {
        Rcpp::stop("Non-finite correlation encountered for SNP pair (" +
                   kept_snps[i]->GetRsid() + ", " + kept_snps[j]->GetRsid() + ").");
      }
      cor_mat(i, j) = cor;
      cor_mat(j, i) = cor;
    }
  }

  Rcpp::StringVector rsid_vec;
  Rcpp::IntegerVector chr_vec;
  Rcpp::IntegerVector bp_vec;
  Rcpp::StringVector a1_vec;
  Rcpp::StringVector a2_vec;
  Rcpp::NumericVector af1pop_vec;

  for (int i = 0; i < m; ++i) {
    Snp* snp = kept_snps[i];
    rsid_vec.push_back(snp->GetRsid());
    chr_vec.push_back(snp->GetChr());
    bp_vec.push_back(snp->GetBp());
    a1_vec.push_back(snp->GetA1());
    a2_vec.push_back(snp->GetA2());
    af1pop_vec.push_back(kept_af1[i]);
  }

  Rcpp::DataFrame snplist = Rcpp::DataFrame::create(
    Rcpp::Named("rsid") = rsid_vec,
    Rcpp::Named("chr") = chr_vec,
    Rcpp::Named("bp") = bp_vec,
    Rcpp::Named("a1") = a1_vec,
    Rcpp::Named("a2") = a2_vec,
    Rcpp::Named("af1pop") = af1pop_vec
  );

  return Rcpp::List::create(
    Rcpp::Named("snplist") = snplist,
    Rcpp::Named("cormat") = cor_mat
  );
}

//' Read superpopulation genotype dosages for one reference window
//'
//' Returns the 0/1/2 dosages (copies of `a1`) of the selected individuals for
//' every reference SNP in the window, with no MAF or missingness filtering.
//' Used to give simulated GWAS and scoring the exact individuals of one part of
//' a sample split (see \code{\link{split_reference_samples}}).
//'
//' @inheritParams computeSuperPopLDWindow
//' @return A list with `snplist` (data frame: rsid, chr, bp, a1, a2),
//'   `genotypes` (integer matrix, individuals x SNPs, `NA` where the reference
//'   has no 0/1/2 call), and `samples` (data frame: population, position --
//'   1-based within the population -- one row per matrix row).
// [[Rcpp::export]]
Rcpp::List readSuperPopGenotypes(
    int chr,
    long long int start_bp,
    long long int end_bp,
    std::string superpopulation_label,
    std::string reference_index_file,
    std::string reference_data_file,
    std::string reference_pop_desc_file,
    Rcpp::Nullable<Rcpp::IntegerVector> keep_bp = R_NilValue,
    Rcpp::Nullable<Rcpp::List> sample_keep = R_NilValue) {

  Arguments args = make_window_args(chr, start_bp, end_bp, superpopulation_label,
                                    reference_index_file, reference_data_file,
                                    reference_pop_desc_file);

  std::map<MapKey, Snp*, LessThanMapKey> snp_map;
  std::vector<Snp*> snp_vec;
  SuperPopLDWindowCleanup cleanup(snp_map, snp_vec);
  std::vector<std::string> selected_pops;
  load_superpop_window(args, snp_map, snp_vec, cleanup, keep_bp, sample_keep, selected_pops);

  // Row labels, from the kept positions (or every position) of each population.
  Rcpp::CharacterVector sample_pop;
  Rcpp::IntegerVector sample_pos;
  const std::vector<std::vector<int> > keep_samples = resolve_sample_keep(args, sample_keep);
  int pop_rank = 0;
  for (int i = 0; i < args.num_pops; ++i) {
    if (!args.pop_flag_vec[i]) continue;
    if (keep_samples.empty()) {
      for (int k = 0; k < args.ref_pop_size_vec[i]; ++k) {
        sample_pop.push_back(args.ref_pop_vec[i]);
        sample_pos.push_back(k + 1);
      }
    } else {
      for (size_t k = 0; k < keep_samples[pop_rank].size(); ++k) {
        sample_pop.push_back(args.ref_pop_vec[i]);
        sample_pos.push_back(keep_samples[pop_rank][k] + 1);
      }
    }
    pop_rank++;
  }

  const int n = sample_pop.size();
  const int m = static_cast<int>(snp_vec.size());
  Rcpp::IntegerMatrix geno(n, m);

  for (int j = 0; j < m; ++j) {
    const std::vector<std::string>& geno_vec = snp_vec[j]->GetGenotypeVec();
    int row = 0;
    for (size_t p = 0; p < geno_vec.size(); ++p) {
      const std::string& g = geno_vec[p];
      for (size_t k = 0; k < g.size(); ++k) {
        if (row >= n) {
          Rcpp::stop("Genotype string longer than population size for SNP " + snp_vec[j]->GetRsid() + ".");
        }
        const char c = g[k];
        geno(row++, j) = (c == '0' || c == '1' || c == '2') ? (c - '0') : NA_INTEGER;
      }
    }
    if (row != n) {
      Rcpp::stop("Genotype string shorter than population size for SNP " + snp_vec[j]->GetRsid() + ".");
    }
  }

  Rcpp::StringVector rsid_vec(m), a1_vec(m), a2_vec(m);
  Rcpp::IntegerVector chr_vec(m), bp_vec(m);
  for (int j = 0; j < m; ++j) {
    rsid_vec[j] = snp_vec[j]->GetRsid();
    chr_vec[j] = snp_vec[j]->GetChr();
    bp_vec[j] = snp_vec[j]->GetBp();
    a1_vec[j] = snp_vec[j]->GetA1();
    a2_vec[j] = snp_vec[j]->GetA2();
  }

  return Rcpp::List::create(
    Rcpp::Named("snplist") = Rcpp::DataFrame::create(
      Rcpp::Named("rsid") = rsid_vec, Rcpp::Named("chr") = chr_vec,
      Rcpp::Named("bp") = bp_vec, Rcpp::Named("a1") = a1_vec,
      Rcpp::Named("a2") = a2_vec, Rcpp::Named("stringsAsFactors") = false),
    Rcpp::Named("genotypes") = geno,
    Rcpp::Named("samples") = Rcpp::DataFrame::create(
      Rcpp::Named("population") = sample_pop, Rcpp::Named("position") = sample_pos,
      Rcpp::Named("stringsAsFactors") = false)
  );
}
