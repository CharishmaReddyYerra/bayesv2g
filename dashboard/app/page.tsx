"use client";

import { useMemo, useState } from "react";
import {
  Activity,
  ArrowUpRight,
  CircleHelp,
  Database,
  Dna,
  FlaskConical,
  Microscope,
} from "lucide-react";

import { Button } from "@/components/ui/button";
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from "@/components/ui/tooltip";
import dashboardData from "@/lib/dashboard-data.json";
import { CredibleVariantTrack } from "@/components/dashboard/credible-variant-track";
import { MethodologySection } from "@/components/dashboard/methodology-section";
import { Sort1AuditSection } from "@/components/dashboard/sort1-audit-section";
import { ValidationSection } from "@/components/dashboard/validation-section";

type ProbabilityMode = "absolute" | "conditional";

type Locus = {
  studyLocusId: string;
  locus_label: string;
  chromosome: number;
  position: number;
  lead_variant: string;
  p_value: string;
  beta: number;
  sample_size: number;
  finemapping_method: string;
  credible_set_size_95: number;
  pip_coverage_95: number;
  candidate_count: number;
};

type Gene = {
  studyLocusId: string;
  locus_label: string;
  geneId: string;
  gene_symbol: string;
  gene_name: string;
  gene_tss: number;
  distance_to_tss_bp: number;
  l2g_score_external: number;
  eQtlColocH4Maximum: number;
  pQtlColocH4Maximum: number;
  sQtlColocH4Maximum: number;
  distanceTssMean: number;
  vepMean: number;
  e2gMean: number;
  conditional_posterior_mean: number;
  conditional_posterior_lower: number;
  conditional_posterior_upper: number;
  posterior_mean: number;
  posterior_lower: number;
  posterior_upper: number;
  unresolved_probability: number;
  candidate_set_is_singleton: boolean;
  probability_scope: string;
  rank: number;
  evidence_coverage: number;
  locus_entropy: number;
};

type DashboardData = {
  release: string;
  loci: Locus[];
  prioritizedGenes: Gene[];
  credibleVariants: Array<{
    studyLocusId: string;
    rsid: string;
    chromosome: number;
    position: number;
    reference_allele: string;
    alternate_allele: string;
    posterior_inclusion_probability: number;
    is_95_credible_set: boolean;
    is_99_credible_set: boolean;
    most_severe_consequence: string;
  }>;
  validation: Array<{
    model: string;
    n_loci: number;
    top1_accuracy: number;
    recall_at_3: number;
    mean_reciprocal_rank: number;
    multiclass_log_loss: number;
    brier_score: number;
    mean_entropy: number;
  }>;
  calibration: Array<{
    bin: string;
    mean_confidence: number;
    observed_accuracy: number;
    n: number;
  }>;
  comparison: {
    top1_difference: number;
    ci_lower_95: number;
    ci_upper_95: number;
    improved_loci: number;
    worsened_loci: number;
    tied_loci: number;
    bootstrap_iterations: number;
  };
  coefficients: Array<{
    feature: string;
    posterior_mean: number;
    posterior_lower: number;
    posterior_upper: number;
    odds_ratio: number;
  }>;
  sort1Colocalization: Array<{
    evidence_type: string;
    dataset_id: string;
    gene_symbol: string;
    shared_snv_count: number;
    pp_h3: number;
    pp_h4: number;
    pp_h4_sensitivity_min: number;
    pp_h4_sensitivity_max: number;
    top_shared_variant: string;
    top_shared_pip_given_h4: number;
    shared_credible_set_size_95: number;
  }>;
  sort1FineMapping: Array<{
    base_pair_location: number;
    effect_allele: string;
    other_allele: string;
    beta: number;
    standard_error: number;
    effect_allele_frequency: number;
    p_value: string;
    rsid: string;
    approximate_pip: number;
    in_95_credible_set: boolean;
    pip_rank: number;
  }>;
  sensitivity: Array<{
    scenario: string;
    top_gene: string;
    top_probability: number;
    reference_gene: string;
    reference_probability: number;
    top_gene_changed: boolean;
    max_probability_change: number;
  }>;
};

const data = dashboardData as unknown as DashboardData;

const probabilityFields = {
  absolute: {
    mean: "posterior_mean",
    lower: "posterior_lower",
    upper: "posterior_upper",
    label: "Absolute probability",
  },
  conditional: {
    mean: "conditional_posterior_mean",
    lower: "conditional_posterior_lower",
    upper: "conditional_posterior_upper",
    label: "Within surfaced candidates",
  },
} as const;

function pct(value: number, digits = 1) {
  return `${(value * 100).toFixed(digits)}%`;
}

function compactPosition(value: number) {
  return `${(value / 1_000_000).toFixed(2)} Mb`;
}

function EvidenceTrack({
  label,
  value,
  note,
}: {
  label: string;
  value: number;
  note: string;
}) {
  return (
    <div className="grid grid-cols-[72px_1fr_42px] items-center gap-3">
      <span className="text-[11px] font-medium text-[var(--ink-muted)]">{label}</span>
      <Tooltip>
        <TooltipTrigger asChild>
          <div className="relative h-[6px] cursor-help overflow-hidden bg-[var(--rule-soft)]">
            <div
              className="absolute inset-y-0 left-0 bg-[var(--cobalt)] transition-[width] duration-300"
              style={{ width: `${Math.max(value * 100, value > 0 ? 2 : 0)}%` }}
            />
          </div>
        </TooltipTrigger>
        <TooltipContent side="top" className="max-w-56 bg-[var(--ink)] text-white">
          {note}
        </TooltipContent>
      </Tooltip>
      <span className="font-mono text-[10px] tabular-nums text-[var(--ink-muted)]">
        {value.toFixed(2)}
      </span>
    </div>
  );
}

function GeneLocusStrip({ locus, genes, selectedGeneId }: { locus: Locus; genes: Gene[]; selectedGeneId: string }) {
  const points = [locus.position, ...genes.map((gene) => gene.gene_tss)];
  const margin = Math.max(15_000, Math.max(...points) - Math.min(...points)) * 0.12;
  const min = Math.min(...points) - margin;
  const max = Math.max(...points) + margin;
  const x = (position: number) => 4 + ((position - min) / (max - min)) * 92;

  return (
    <div className="mt-6 border-y border-[var(--rule)] py-4">
      <div className="mb-3 flex items-center justify-between text-[10px] uppercase tracking-[0.12em] text-[var(--ink-faint)]">
        <span>Local gene map</span>
        <span className="font-mono normal-case tracking-normal">chr{locus.chromosome}:{compactPosition(locus.position)}</span>
      </div>
      <div className="relative h-14" aria-label={`Gene positions around ${locus.lead_variant}`}>
        <div className="absolute left-[4%] right-[4%] top-6 h-px bg-[var(--rule-strong)]" />
        <div className="absolute top-[18px] h-4 w-px bg-[var(--coral)]" style={{ left: `${x(locus.position)}%` }}>
          <span className="absolute -left-8 -top-5 w-16 text-center font-mono text-[9px] text-[var(--coral-deep)]">
            {locus.lead_variant}
          </span>
        </div>
        {genes.map((gene, index) => (
          <div
            key={gene.geneId}
            className="absolute top-5 -translate-x-1/2"
            style={{ left: `${x(gene.gene_tss)}%` }}
          >
            <div
              className={`mx-auto h-3 w-3 rotate-45 border ${
                gene.geneId === selectedGeneId
                  ? "border-[var(--cobalt)] bg-[var(--cobalt)]"
                  : "border-[var(--ink-faint)] bg-[var(--paper)]"
              }`}
            />
            <span className={`mt-2 block -translate-x-[calc(50%-6px)] whitespace-nowrap text-[9px] ${index % 2 ? "pt-3" : ""}`}>
              {gene.gene_symbol}
            </span>
          </div>
        ))}
      </div>
    </div>
  );
}

export default function Home() {
  const initialLocus = data.loci.find((locus) => locus.locus_label.includes("SORT1")) ?? data.loci[0];
  const initialGene = data.prioritizedGenes.find(
    (gene) => gene.studyLocusId === initialLocus.studyLocusId && gene.rank === 1,
  )!;

  const [selectedLocusId, setSelectedLocusId] = useState(initialLocus.studyLocusId);
  const [selectedGeneId, setSelectedGeneId] = useState(initialGene.geneId);
  const [probabilityMode, setProbabilityMode] = useState<ProbabilityMode>("absolute");

  const selectedLocus = data.loci.find((locus) => locus.studyLocusId === selectedLocusId)!;
  const locusGenes = useMemo(
    () =>
      data.prioritizedGenes
        .filter((gene) => gene.studyLocusId === selectedLocusId)
        .sort((a, b) => a.rank - b.rank),
    [selectedLocusId],
  );
  const selectedGene = locusGenes.find((gene) => gene.geneId === selectedGeneId) ?? locusGenes[0];
  const topGeneByLocus = useMemo(
    () => new Map(data.prioritizedGenes.filter((gene) => gene.rank === 1).map((gene) => [gene.studyLocusId, gene])),
    [],
  );
  const validation = data.validation.find((row) => row.model === "BayesV2G")!;
  const field = probabilityFields[probabilityMode];
  const selectedVariants = data.credibleVariants.filter((variant) => variant.studyLocusId === selectedLocusId);
  const sort1LocusId = data.loci.find((locus) => locus.locus_label.includes("SORT1"))!.studyLocusId;
  const sort1Genes = data.prioritizedGenes.filter((gene) => gene.studyLocusId === sort1LocusId);

  const selectLocus = (locusId: string) => {
    const topGene = topGeneByLocus.get(locusId)!;
    setSelectedLocusId(locusId);
    setSelectedGeneId(topGene.geneId);
  };

  const selectedProbability = selectedGene[field.mean];
  const interpretation = selectedGene.candidate_set_is_singleton
    ? "Singleton candidate set: the within-set probability is structurally 100%, so use the absolute probability and unresolved mass for interpretation."
    : selectedGene.evidence_coverage < 0.5
      ? "Sparse evidence profile. A zero can mean no signal, unavailable tissue, or low power—not evidence against the gene."
      : `${selectedGene.gene_symbol} is supported by ${Math.round(selectedGene.evidence_coverage * 5)} of 5 modeled evidence families in this snapshot.`;

  return (
    <TooltipProvider>
      <main className="min-h-screen bg-[var(--canvas)] text-[var(--ink)]">
        <header className="sticky top-0 z-30 border-b border-[var(--rule)] bg-[color:var(--paper)/0.96] backdrop-blur-sm">
          <div className="mx-auto flex h-16 max-w-[1580px] items-center justify-between px-4 sm:px-6 lg:px-8">
            <div className="flex items-center gap-3">
              <div className="grid size-9 place-items-center border border-[var(--ink)] bg-[var(--ink)] text-white">
                <Dna className="size-[18px]" strokeWidth={1.7} />
              </div>
              <div>
                <p className="text-[10px] font-semibold uppercase tracking-[0.18em] text-[var(--cobalt)]">BayesV2G</p>
                <h1 className="text-[15px] font-semibold tracking-[-0.01em]">Evidence Workbench</h1>
              </div>
            </div>
            <div className="hidden items-center gap-7 md:flex">
              <a href="#explorer" className="text-xs font-medium text-[var(--ink)] hover:text-[var(--cobalt)]">Locus explorer</a>
              <a href="#validation" className="text-xs text-[var(--ink-faint)] hover:text-[var(--cobalt)]">Validation</a>
              <a href="#sort1-audit" className="text-xs text-[var(--ink-faint)] hover:text-[var(--cobalt)]">SORT1 audit</a>
              <a href="#methods" className="text-xs text-[var(--ink-faint)] hover:text-[var(--cobalt)]">Methods</a>
              <span className="h-4 w-px bg-[var(--rule)]" />
              <span className="flex items-center gap-2 font-mono text-[10px] text-[var(--ink-muted)]">
                <span className="size-1.5 bg-[var(--teal)]" /> Real data · {data.release}
              </span>
            </div>
          </div>
        </header>

        <section id="explorer" className="mx-auto max-w-[1580px] scroll-mt-20 px-4 py-5 sm:px-6 lg:px-8">
          <div className="mb-5 flex flex-col justify-between gap-4 border-b border-[var(--rule)] pb-5 lg:flex-row lg:items-end">
            <div>
              <div className="mb-2 flex items-center gap-2 text-[10px] font-semibold uppercase tracking-[0.16em] text-[var(--ink-faint)]">
                <FlaskConical className="size-3.5" /> LDL cholesterol · ten-locus demonstration
              </div>
              <h2 className="font-serif text-2xl leading-tight tracking-[-0.025em] sm:text-[31px]">
                Resolve a locus into a ranked, auditable gene hypothesis.
              </h2>
            </div>
            <p className="max-w-xl text-xs leading-5 text-[var(--ink-muted)]">
              Bayesian integration of fine-mapping, molecular QTL colocalization, functional consequence,
              enhancer-to-gene evidence, and distance. Probabilities retain unresolved mass.
            </p>
          </div>

          <div className="grid gap-4 xl:grid-cols-[230px_minmax(0,1fr)_330px]">
            <aside className="border border-[var(--rule)] bg-[var(--paper)]">
              <div className="flex h-11 items-center justify-between border-b border-[var(--rule)] px-3.5">
                <span className="text-[10px] font-semibold uppercase tracking-[0.14em] text-[var(--ink-muted)]">LDL loci</span>
                <span className="font-mono text-[10px] text-[var(--ink-faint)]">10</span>
              </div>
              <nav aria-label="LDL locus selection" className="max-h-[300px] overflow-y-auto xl:max-h-[650px]">
                {data.loci.map((locus, index) => {
                  const topGene = topGeneByLocus.get(locus.studyLocusId)!;
                  const active = locus.studyLocusId === selectedLocusId;
                  return (
                    <button
                      key={locus.studyLocusId}
                      type="button"
                      onClick={() => selectLocus(locus.studyLocusId)}
                      className={`group grid w-full grid-cols-[24px_1fr] gap-2 border-b border-[var(--rule-soft)] px-3 py-3 text-left transition-colors last:border-b-0 ${
                        active ? "bg-[var(--cobalt-wash)]" : "hover:bg-[var(--canvas)]"
                      }`}
                    >
                      <span className={`font-mono text-[10px] ${active ? "text-[var(--cobalt)]" : "text-[var(--ink-faint)]"}`}>
                        {String(index + 1).padStart(2, "0")}
                      </span>
                      <span className="min-w-0">
                        <span className="block truncate text-xs font-semibold">{locus.locus_label.replace(" region", "")}</span>
                        <span className="mt-1 flex items-center justify-between gap-2 text-[10px] text-[var(--ink-muted)]">
                          <span className="font-mono">chr{locus.chromosome}:{compactPosition(locus.position)}</span>
                          <span className={active ? "text-[var(--cobalt)]" : ""}>{topGene.gene_symbol}</span>
                        </span>
                      </span>
                    </button>
                  );
                })}
              </nav>
            </aside>

            <section className="min-w-0 border border-[var(--rule)] bg-[var(--paper)]">
              <div className="border-b border-[var(--rule)] px-4 py-4 sm:px-5">
                <div className="flex flex-col justify-between gap-4 md:flex-row md:items-start">
                  <div>
                    <div className="flex flex-wrap items-center gap-2">
                      <h3 className="text-lg font-semibold tracking-[-0.02em]">{selectedLocus.locus_label}</h3>
                      <span className="border border-[var(--rule-strong)] px-1.5 py-0.5 font-mono text-[9px] text-[var(--ink-muted)]">
                        {selectedLocus.lead_variant}
                      </span>
                    </div>
                    <p className="mt-1 font-mono text-[10px] text-[var(--ink-faint)]">
                      chr{selectedLocus.chromosome}:{selectedLocus.position.toLocaleString()} · β {selectedLocus.beta.toFixed(3)} · n={selectedLocus.sample_size.toLocaleString()}
                    </p>
                  </div>
                  <div className="inline-flex w-fit border border-[var(--rule-strong)] p-0.5" aria-label="Probability display mode">
                    {(["absolute", "conditional"] as ProbabilityMode[]).map((mode) => (
                      <Button
                        key={mode}
                        size="xs"
                        variant="ghost"
                        onClick={() => setProbabilityMode(mode)}
                        className={`h-7 rounded-none px-2.5 text-[10px] ${
                          probabilityMode === mode
                            ? "bg-[var(--ink)] text-white hover:bg-[var(--ink)] hover:text-white"
                            : "text-[var(--ink-muted)]"
                        }`}
                      >
                        {mode === "absolute" ? "Absolute" : "Within set"}
                      </Button>
                    ))}
                  </div>
                </div>
                <div className="mt-4 grid grid-cols-3 gap-px border border-[var(--rule)] bg-[var(--rule)]">
                  <div className="bg-[var(--canvas)] px-3 py-2.5">
                    <p className="text-[9px] uppercase tracking-[0.12em] text-[var(--ink-faint)]">95% credible set</p>
                    <p className="mt-1 font-mono text-sm font-semibold">{selectedLocus.credible_set_size_95} variants</p>
                  </div>
                  <div className="bg-[var(--canvas)] px-3 py-2.5">
                    <p className="text-[9px] uppercase tracking-[0.12em] text-[var(--ink-faint)]">PIP captured</p>
                    <p className="mt-1 font-mono text-sm font-semibold">{pct(selectedLocus.pip_coverage_95, 2)}</p>
                  </div>
                  <div className="bg-[var(--canvas)] px-3 py-2.5">
                    <p className="text-[9px] uppercase tracking-[0.12em] text-[var(--ink-faint)]">Candidates</p>
                    <p className="mt-1 font-mono text-sm font-semibold">{selectedLocus.candidate_count} surfaced</p>
                  </div>
                </div>
              </div>

              <div className="px-4 py-5 sm:px-5">
                <div className="mb-4 flex items-end justify-between gap-4">
                  <div>
                    <p className="text-[10px] font-semibold uppercase tracking-[0.14em] text-[var(--ink-muted)]">Ranked candidates</p>
                    <p className="mt-1 text-[10px] text-[var(--ink-faint)]">{field.label} · posterior mean with 90% interval</p>
                  </div>
                  {probabilityMode === "absolute" && (
                    <span className="text-[10px] text-[var(--ink-faint)]">
                      unresolved <strong className="font-mono font-medium text-[var(--coral-deep)]">{pct(locusGenes[0].unresolved_probability)}</strong>
                    </span>
                  )}
                </div>

                <div className="space-y-2">
                  {locusGenes.map((gene) => {
                    const active = gene.geneId === selectedGene.geneId;
                    const mean = gene[field.mean];
                    const lower = gene[field.lower];
                    const upper = gene[field.upper];
                    return (
                      <button
                        key={gene.geneId}
                        type="button"
                        onClick={() => setSelectedGeneId(gene.geneId)}
                        className={`grid w-full grid-cols-[minmax(86px,0.32fr)_minmax(150px,1fr)_56px] items-center gap-3 border px-3 py-2.5 text-left transition-colors ${
                          active
                            ? "border-[var(--cobalt)] bg-[var(--cobalt-wash)]"
                            : "border-transparent hover:border-[var(--rule)] hover:bg-[var(--canvas)]"
                        }`}
                      >
                        <span>
                          <span className="block text-xs font-semibold">{gene.gene_symbol}</span>
                          <span className="mt-0.5 block font-mono text-[9px] text-[var(--ink-faint)]">
                            {Math.round(gene.distance_to_tss_bp / 1000)} kb to TSS
                          </span>
                        </span>
                        <span className="relative block h-2 bg-[var(--rule-soft)]">
                          <span
                            className={`absolute inset-y-0 left-0 ${active ? "bg-[var(--cobalt)]" : "bg-[var(--slate)]"}`}
                            style={{ width: `${mean * 100}%` }}
                          />
                          <span
                            className="absolute top-1/2 h-[14px] -translate-y-1/2 border-x border-[var(--ink-muted)]"
                            style={{ left: `${lower * 100}%`, width: `${Math.max((upper - lower) * 100, 0.4)}%` }}
                          />
                        </span>
                        <span className="text-right font-mono text-xs font-semibold tabular-nums">{pct(mean)}</span>
                      </button>
                    );
                  })}
                </div>

                <GeneLocusStrip locus={selectedLocus} genes={locusGenes} selectedGeneId={selectedGene.geneId} />
                <CredibleVariantTrack
                  key={selectedLocusId}
                  variants={selectedVariants}
                  leadVariant={selectedLocus.lead_variant}
                />
              </div>
            </section>

            <aside className="border border-[var(--rule)] bg-[var(--paper)]">
              <div className="border-b border-[var(--rule)] px-4 py-4">
                <div className="flex items-start justify-between gap-3">
                  <div>
                    <p className="text-[10px] font-semibold uppercase tracking-[0.14em] text-[var(--cobalt)]">Evidence inspector</p>
                    <h3 className="mt-1 text-xl font-semibold tracking-[-0.025em]">{selectedGene.gene_symbol}</h3>
                  </div>
                  <span className="grid size-8 place-items-center border border-[var(--rule)] font-mono text-[10px]">#{selectedGene.rank}</span>
                </div>
                <p className="mt-2 text-xs leading-5 text-[var(--ink-muted)]">{selectedGene.gene_name}</p>
              </div>

              <div className="border-b border-[var(--rule)] px-4 py-4">
                <div className="flex items-end justify-between">
                  <div>
                    <p className="text-[9px] uppercase tracking-[0.12em] text-[var(--ink-faint)]">{field.label}</p>
                    <p className="mt-1 font-serif text-3xl tracking-[-0.03em]">{pct(selectedProbability)}</p>
                  </div>
                  <div className="text-right font-mono text-[9px] leading-4 text-[var(--ink-faint)]">
                    <span>{pct(selectedGene[field.lower])}</span>
                    <span className="mx-1">—</span>
                    <span>{pct(selectedGene[field.upper])}</span>
                    <br />90% interval
                  </div>
                </div>
              </div>

              <div className="space-y-4 px-4 py-5">
                <div className="flex items-center justify-between">
                  <p className="text-[10px] font-semibold uppercase tracking-[0.14em] text-[var(--ink-muted)]">Modeled evidence</p>
                  <Tooltip>
                    <TooltipTrigger asChild>
                      <button type="button" aria-label="Explain evidence scale" className="text-[var(--ink-faint)] hover:text-[var(--ink)]">
                        <CircleHelp className="size-3.5" />
                      </button>
                    </TooltipTrigger>
                    <TooltipContent className="max-w-60 bg-[var(--ink)] text-white">
                      Tracks are normalized feature values, not independent probabilities.
                    </TooltipContent>
                  </Tooltip>
                </div>
                <EvidenceTrack label="Proximity" value={selectedGene.distanceTssMean} note="Distance-to-TSS feature; higher is nearer on the model scale." />
                <EvidenceTrack label="eQTL H4" value={selectedGene.eQtlColocH4Maximum} note="Maximum posterior probability of a shared eQTL and GWAS causal variant." />
                <EvidenceTrack label="pQTL H4" value={selectedGene.pQtlColocH4Maximum} note="Maximum posterior probability of pQTL colocalization." />
                <EvidenceTrack label="sQTL H4" value={selectedGene.sQtlColocH4Maximum} note="Maximum posterior probability of sQTL colocalization." />
                <EvidenceTrack label="VEP" value={selectedGene.vepMean} note="Aggregated functional consequence evidence from VEP annotations." />
                <EvidenceTrack label="E2G" value={selectedGene.e2gMean} note="Enhancer-to-gene linking evidence in the locus." />

                <div className="mt-5 border-t border-[var(--rule)] pt-4">
                  <div className="flex items-center justify-between text-[10px]">
                    <span className="text-[var(--ink-muted)]">External Open Targets L2G</span>
                    <span className="font-mono font-semibold">{selectedGene.l2g_score_external.toFixed(3)}</span>
                  </div>
                  <p className="mt-2 text-[9px] leading-4 text-[var(--ink-faint)]">
                    Used to surface candidates (score ≥ 0.05), not as a fitted BayesV2G feature.
                  </p>
                </div>
              </div>

              <div className="m-4 border-l-2 border-[var(--coral)] bg-[var(--coral-wash)] px-3 py-3">
                <p className="text-[9px] font-semibold uppercase tracking-[0.13em] text-[var(--coral-deep)]">Interpretation note</p>
                <p className="mt-1.5 text-[10px] leading-[1.55] text-[var(--ink-muted)]">{interpretation}</p>
              </div>
            </aside>
          </div>

          <div className="mt-4 grid border border-[var(--rule)] bg-[var(--paper)] md:grid-cols-[1.1fr_repeat(3,1fr)]">
            <div className="flex items-center gap-3 border-b border-[var(--rule)] px-4 py-3 md:border-b-0 md:border-r">
              <Activity className="size-4 text-[var(--cobalt)]" />
              <div>
                <p className="text-[10px] font-semibold uppercase tracking-[0.12em]">Held-out validation</p>
                <p className="text-[9px] text-[var(--ink-faint)]">1,048 independent Open Targets loci</p>
              </div>
            </div>
            <div className="border-b border-[var(--rule)] px-4 py-3 md:border-b-0 md:border-r">
              <p className="text-[9px] uppercase tracking-[0.11em] text-[var(--ink-faint)]">Top-1 accuracy</p>
              <p className="mt-1 font-mono text-lg font-semibold">{pct(validation.top1_accuracy, 2)}</p>
            </div>
            <div className="border-b border-[var(--rule)] px-4 py-3 md:border-b-0 md:border-r">
              <p className="text-[9px] uppercase tracking-[0.11em] text-[var(--ink-faint)]">Recall @ 3</p>
              <p className="mt-1 font-mono text-lg font-semibold">{pct(validation.recall_at_3, 2)}</p>
            </div>
            <div className="px-4 py-3">
              <p className="text-[9px] uppercase tracking-[0.11em] text-[var(--ink-faint)]">Gain vs distance</p>
              <p className="mt-1 flex items-center gap-1 font-mono text-lg font-semibold text-[var(--teal-deep)]">
                +{pct(data.comparison.top1_difference, 2)} <ArrowUpRight className="size-3.5" />
              </p>
            </div>
          </div>

          <div className="mt-3 flex flex-col justify-between gap-2 text-[9px] leading-4 text-[var(--ink-faint)] sm:flex-row">
            <span className="flex items-center gap-1.5"><Database className="size-3" /> Real public genetics data · reproducible project outputs</span>
            <span className="flex items-center gap-1.5"><Microscope className="size-3" /> Research prioritization only—not clinical guidance</span>
          </div>

          <div className="mt-12 space-y-12">
            <ValidationSection
              validation={data.validation}
              calibration={data.calibration}
              comparison={data.comparison}
              coefficients={data.coefficients}
            />
            <Sort1AuditSection
              colocalization={data.sort1Colocalization}
              fineMapping={data.sort1FineMapping}
              genes={sort1Genes}
              sensitivity={data.sensitivity}
            />
            <MethodologySection />
          </div>
        </section>
      </main>
    </TooltipProvider>
  );
}
