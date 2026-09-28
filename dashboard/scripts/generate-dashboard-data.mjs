import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const resultsDir = path.resolve(here, "../../results/ldl_demo");
const outputPath = path.resolve(here, "../lib/dashboard-data.json");

function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;

  for (let index = 0; index < text.length; index += 1) {
    const char = text[index];
    const next = text[index + 1];

    if (char === '"' && quoted && next === '"') {
      field += '"';
      index += 1;
    } else if (char === '"') {
      quoted = !quoted;
    } else if (char === "," && !quoted) {
      row.push(field);
      field = "";
    } else if ((char === "\n" || char === "\r") && !quoted) {
      if (char === "\r" && next === "\n") index += 1;
      row.push(field);
      if (row.some((value) => value.length > 0)) rows.push(row);
      row = [];
      field = "";
    } else {
      field += char;
    }
  }

  if (field.length > 0 || row.length > 0) {
    row.push(field);
    rows.push(row);
  }

  const [headers, ...values] = rows;
  return values.map((cells) =>
    Object.fromEntries(headers.map((header, index) => [header, cells[index] ?? ""])),
  );
}

function readCsv(filename) {
  return parseCsv(fs.readFileSync(path.join(resultsDir, filename), "utf8"));
}

function numeric(value) {
  if (value === "" || value == null) return null;
  const parsed = Number(value);
  return Number.isFinite(parsed) ? parsed : value;
}

function coerce(value) {
  if (value === "True") return true;
  if (value === "False") return false;
  if (/^-?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?$/i.test(value ?? "")) {
    return numeric(value);
  }
  return value;
}

function pick(row, fields) {
  return Object.fromEntries(
    fields.map((field) => [
      field,
      field === "p_value" ? row[field] : coerce(row[field]),
    ]),
  );
}

const prioritizedFields = [
  "studyLocusId",
  "locus_label",
  "lead_variant",
  "chromosome",
  "position",
  "geneId",
  "gene_symbol",
  "gene_name",
  "gene_biotype",
  "gene_tss",
  "distance_to_tss_bp",
  "l2g_score_external",
  "eQtlColocH4Maximum",
  "pQtlColocH4Maximum",
  "sQtlColocH4Maximum",
  "distanceTssMean",
  "vepMean",
  "e2gMean",
  "geneCount500kb",
  "proteinGeneCount500kb",
  "conditional_posterior_mean",
  "conditional_posterior_lower",
  "conditional_posterior_upper",
  "posterior_mean",
  "posterior_lower",
  "posterior_upper",
  "unresolved_probability",
  "surfaced_candidate_count",
  "candidate_set_is_singleton",
  "probability_scope",
  "rank",
  "evidence_coverage",
  "locus_entropy",
];

const validationFields = [
  "model",
  "n_loci",
  "top1_accuracy",
  "recall_at_3",
  "mean_reciprocal_rank",
  "multiclass_log_loss",
  "brier_score",
  "mean_entropy",
];

const dashboard = {
  generatedFrom: "BayesV2G LDL demo outputs",
  release: "Open Targets 26.09 / GTEx v8-v10",
  loci: readCsv("loci.csv").map((row) =>
    pick(row, [
      "studyLocusId",
      "locus_label",
      "study_id",
      "trait",
      "chromosome",
      "position",
      "lead_variant",
      "p_value",
      "beta",
      "sample_size",
      "finemapping_method",
      "credible_set_confidence",
      "credible_set_size_95",
      "pip_coverage_95",
      "candidate_count",
      "publication_title",
      "pubmed_id",
    ]),
  ),
  prioritizedGenes: readCsv("prioritized_genes.csv").map((row) => pick(row, prioritizedFields)),
  credibleVariants: readCsv("credible_set_variants.csv").map((row) =>
    pick(row, [
      "studyLocusId",
      "rsid",
      "chromosome",
      "position",
      "reference_allele",
      "alternate_allele",
      "posterior_inclusion_probability",
      "is_95_credible_set",
      "is_99_credible_set",
      "log_bayes_factor",
      "most_severe_consequence",
    ]),
  ),
  validation: readCsv("validation_metrics.csv").map((row) => pick(row, validationFields)),
  calibration: readCsv("calibration.csv").map((row) =>
    pick(row, ["bin", "mean_confidence", "observed_accuracy", "n"]),
  ),
  comparison: readCsv("paired_top1_comparison.csv").map((row) =>
    pick(row, [
      "first_model",
      "second_model",
      "n_loci",
      "top1_difference",
      "ci_lower_95",
      "ci_upper_95",
      "improved_loci",
      "worsened_loci",
      "tied_loci",
      "bootstrap_iterations",
    ]),
  )[0],
  coefficients: readCsv("posterior_coefficients.csv").map((row) =>
    pick(row, ["feature", "posterior_mean", "posterior_lower", "posterior_upper", "odds_ratio"]),
  ),
  sort1Colocalization: readCsv("sort1_colocalization.csv").map((row) =>
    pick(row, [
      "evidence_type",
      "dataset_id",
      "gene_symbol",
      "shared_snv_count",
      "pp_h3",
      "pp_h4",
      "pp_h4_sensitivity_min",
      "pp_h4_sensitivity_max",
      "top_shared_variant",
      "top_shared_pip_given_h4",
      "shared_credible_set_size_95",
    ]),
  ),
  sort1FineMapping: readCsv("sort1_finemapping_top25.csv").map((row) =>
    pick(row, [
      "base_pair_location",
      "effect_allele",
      "other_allele",
      "beta",
      "standard_error",
      "effect_allele_frequency",
      "p_value",
      "rsid",
      "approximate_pip",
      "in_95_credible_set",
      "pip_rank",
    ]),
  ),
  sensitivity: readCsv("sensitivity.csv")
    .filter((row) => row.studyLocusId === "5eb8878b948ed248f0832a6a1b27d324")
    .map((row) =>
      pick(row, [
        "scenario",
        "scenario_type",
        "top_gene",
        "top_probability",
        "reference_gene",
        "reference_probability",
        "top_gene_changed",
        "max_probability_change",
      ]),
    ),
};

fs.writeFileSync(outputPath, `${JSON.stringify(dashboard, null, 2)}\n`);
console.log(`Wrote ${outputPath}`);
