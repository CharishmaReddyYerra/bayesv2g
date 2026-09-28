"use client";

import { useMemo, useState } from "react";
import { Braces, CheckCircle2, FlaskConical, TriangleAlert } from "lucide-react";

type Coloc = {
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
};

type FineMapVariant = {
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
};

type Gene = {
  gene_symbol: string;
  gene_name: string;
  rank: number;
  posterior_mean: number;
  conditional_posterior_mean: number;
  l2g_score_external: number;
  eQtlColocH4Maximum: number;
  pQtlColocH4Maximum: number;
  sQtlColocH4Maximum: number;
  vepMean: number;
  e2gMean: number;
};

type Sensitivity = {
  scenario: string;
  top_gene: string;
  top_probability: number;
  reference_gene: string;
  reference_probability: number;
  top_gene_changed: boolean;
  max_probability_change: number;
};

function pct(value: number, digits = 2) {
  return `${(value * 100).toFixed(digits)}%`;
}

function EvidenceCell({ value }: { value: number }) {
  return (
    <span className="relative block h-5 min-w-12 border border-[var(--rule-soft)] bg-[var(--canvas)]" title={value.toFixed(4)}>
      <span className="absolute inset-y-0 left-0 bg-[var(--cobalt-wash)]" style={{ width: `${Math.max(value * 100, value > 0 ? 3 : 0)}%` }} />
      <span className="relative grid h-full place-items-center font-mono text-[8px] text-[var(--ink-muted)]">{value.toFixed(2)}</span>
    </span>
  );
}

export function Sort1AuditSection({
  colocalization,
  fineMapping,
  genes,
  sensitivity,
}: {
  colocalization: Coloc[];
  fineMapping: FineMapVariant[];
  genes: Gene[];
  sensitivity: Sensitivity[];
}) {
  const [selectedColocKey, setSelectedColocKey] = useState("liver_eQTL-SORT1");
  const selectedColoc = colocalization.find(
    (row) => `${row.evidence_type}-${row.gene_symbol}` === selectedColocKey,
  ) ?? colocalization[0];
  const rankedVariants = useMemo(
    () => [...fineMapping].sort((a, b) => a.pip_rank - b.pip_rank),
    [fineMapping],
  );
  const minPosition = Math.min(...fineMapping.map((row) => row.base_pair_location));
  const maxPosition = Math.max(...fineMapping.map((row) => row.base_pair_location));
  const plotX = (position: number) => 18 + ((position - minPosition) / (maxPosition - minPosition)) * 264;
  const maxSensitivityShift = Math.max(...sensitivity.map((row) => row.max_probability_change));
  const allStable = sensitivity.every((row) => !row.top_gene_changed);

  return (
    <section id="sort1-audit" className="scroll-mt-20 border-t border-[var(--rule-strong)] pt-10">
      <div className="mb-6 grid gap-4 lg:grid-cols-[1fr_0.72fr] lg:items-end">
        <div>
          <p className="mb-2 flex items-center gap-2 text-[10px] font-semibold uppercase tracking-[0.16em] text-[var(--coral-deep)]">
            <FlaskConical className="size-3.5" /> Deep evidence audit · chr1p13
          </p>
          <h2 className="font-serif text-3xl tracking-[-0.03em]">The SORT1 locus, left intentionally complicated.</h2>
        </div>
        <p className="text-xs leading-5 text-[var(--ink-muted)]">
          A real regional GWAS, independent approximate-Bayes-factor fine-mapping, liver eQTL, and plasma pQTL data expose why an auditable ranking matters.
        </p>
      </div>

      <div className="grid gap-4 xl:grid-cols-[1.04fr_0.96fr]">
        <div className="border border-[var(--rule)] bg-[var(--paper)]">
          <div className="flex items-start justify-between gap-4 border-b border-[var(--rule)] px-4 py-3">
            <div>
              <p className="text-[10px] font-semibold uppercase tracking-[0.13em]">Independent fine-mapping audit</p>
              <p className="mt-1 text-[9px] text-[var(--ink-faint)]">GCST90018961 · 25 displayed regional variants</p>
            </div>
            <span className="font-mono text-[9px] text-[var(--ink-muted)]">lead {rankedVariants[0].rsid}</span>
          </div>
          <div className="p-4">
            <div className="grid grid-cols-[1fr_auto] gap-5">
              <svg viewBox="0 0 300 150" className="h-[150px] w-full" role="img" aria-label="SORT1 regional fine-mapping PIP plot">
                <line x1="18" y1="126" x2="282" y2="126" stroke="var(--rule-strong)" />
                <line x1="18" y1="25" x2="282" y2="25" stroke="var(--rule-soft)" strokeDasharray="3 4" />
                <line x1="18" y1="75" x2="282" y2="75" stroke="var(--rule-soft)" strokeDasharray="3 4" />
                {fineMapping.map((variant) => {
                  const height = Math.max(Math.sqrt(variant.approximate_pip) * 95, 2);
                  const active = variant.pip_rank === 1;
                  return (
                    <g key={`${variant.rsid}-${variant.base_pair_location}`}>
                      <line
                        x1={plotX(variant.base_pair_location)}
                        y1="126"
                        x2={plotX(variant.base_pair_location)}
                        y2={126 - height}
                        stroke={active ? "var(--coral)" : variant.in_95_credible_set ? "var(--cobalt)" : "var(--slate)"}
                        strokeWidth={active ? 3 : 1.5}
                      />
                      <circle cx={plotX(variant.base_pair_location)} cy={126 - height} r={active ? 3.5 : 2} fill={active ? "var(--coral)" : "var(--slate)"}>
                        <title>{`${variant.rsid}: PIP ${variant.approximate_pip.toPrecision(4)}`}</title>
                      </circle>
                    </g>
                  );
                })}
                <text x="18" y="143" fill="var(--ink-faint)" fontSize="7" fontFamily="monospace">{minPosition.toLocaleString()}</text>
                <text x="282" y="143" textAnchor="end" fill="var(--ink-faint)" fontSize="7" fontFamily="monospace">{maxPosition.toLocaleString()}</text>
                <text x="8" y="80" transform="rotate(-90 8 80)" textAnchor="middle" fill="var(--ink-faint)" fontSize="7">√PIP scale</text>
              </svg>
              <div className="hidden w-28 border-l border-[var(--rule)] pl-3 sm:block">
                <p className="text-[8px] uppercase tracking-[0.1em] text-[var(--ink-faint)]">Top PIP</p>
                <p className="mt-1 font-serif text-2xl">{pct(rankedVariants[0].approximate_pip)}</p>
                <p className="mt-1 font-mono text-[9px] text-[var(--coral-deep)]">{rankedVariants[0].rsid}</p>
                <p className="mt-4 text-[8px] uppercase tracking-[0.1em] text-[var(--ink-faint)]">95% set</p>
                <p className="mt-1 font-mono text-sm font-semibold">1 variant</p>
              </div>
            </div>
            <div className="mt-4 overflow-x-auto">
              <table className="w-full min-w-[520px] border-collapse text-left">
                <thead className="border-y border-[var(--rule)] text-[8px] uppercase tracking-[0.1em] text-[var(--ink-faint)]">
                  <tr><th className="py-2 pr-3 font-medium">Rank</th><th className="py-2 pr-3 font-medium">Variant</th><th className="py-2 pr-3 font-medium">Alleles</th><th className="py-2 pr-3 font-medium">β ± SE</th><th className="py-2 pr-3 font-medium">P value</th><th className="py-2 text-right font-medium">Approx. PIP</th></tr>
                </thead>
                <tbody>
                  {rankedVariants.slice(0, 4).map((variant) => (
                    <tr key={variant.rsid} className="border-b border-[var(--rule-soft)] text-[9px]">
                      <td className="py-2 pr-3 font-mono text-[var(--ink-faint)]">{String(variant.pip_rank).padStart(2, "0")}</td>
                      <td className="py-2 pr-3 font-mono font-semibold">{variant.rsid}</td>
                      <td className="py-2 pr-3 font-mono text-[var(--ink-muted)]">{variant.other_allele}›{variant.effect_allele}</td>
                      <td className="py-2 pr-3 font-mono text-[var(--ink-muted)]">{variant.beta.toFixed(4)} ± {variant.standard_error.toFixed(4)}</td>
                      <td className="py-2 pr-3 font-mono text-[var(--ink-muted)]">{variant.p_value}</td>
                      <td className="py-2 text-right font-mono font-semibold">{pct(variant.approximate_pip, variant.approximate_pip > 0.01 ? 2 : 4)}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        </div>

        <div className="border border-[var(--rule)] bg-[var(--paper)]">
          <div className="border-b border-[var(--rule)] px-4 py-3">
            <p className="text-[10px] font-semibold uppercase tracking-[0.13em]">Colocalization ledger</p>
            <p className="mt-1 text-[9px] text-[var(--ink-faint)]">Select a row to inspect the H3/H4 evidence balance</p>
          </div>
          <div className="grid lg:grid-cols-[1fr_190px]">
            <div className="border-b border-[var(--rule)] lg:border-b-0 lg:border-r">
              {colocalization.map((row) => {
                const key = `${row.evidence_type}-${row.gene_symbol}`;
                const active = key === `${selectedColoc.evidence_type}-${selectedColoc.gene_symbol}`;
                return (
                  <button
                    key={key}
                    type="button"
                    onClick={() => setSelectedColocKey(key)}
                    className={`grid w-full grid-cols-[1fr_72px] items-center gap-3 border-b border-[var(--rule-soft)] px-4 py-3 text-left last:border-b-0 ${active ? "bg-[var(--cobalt-wash)]" : "hover:bg-[var(--canvas)]"}`}
                  >
                    <span>
                      <span className="flex items-center gap-2 text-[10px] font-semibold">
                        {row.gene_symbol}
                        <span className={`px-1 py-0.5 font-mono text-[7px] ${row.evidence_type === "liver_eQTL" ? "bg-[#e4efe9] text-[var(--teal-deep)]" : "bg-[var(--coral-wash)] text-[var(--coral-deep)]"}`}>
                          {row.evidence_type === "liver_eQTL" ? "LIVER eQTL" : "PLASMA pQTL"}
                        </span>
                      </span>
                      <span className="mt-1 block font-mono text-[8px] text-[var(--ink-faint)]">{row.dataset_id} · {row.shared_snv_count.toLocaleString()} shared SNVs</span>
                    </span>
                    <span className="text-right font-mono text-[10px] font-semibold">{pct(row.pp_h4, row.pp_h4 > 0.1 ? 3 : 2)}</span>
                  </button>
                );
              })}
            </div>
            <div className="px-4 py-4">
              <p className="text-[8px] uppercase tracking-[0.1em] text-[var(--ink-faint)]">PP(H4): shared causal</p>
              <p className={`mt-1 font-serif text-3xl ${selectedColoc.pp_h4 > 0.8 ? "text-[var(--teal-deep)]" : "text-[var(--coral-deep)]"}`}>
                {pct(selectedColoc.pp_h4, selectedColoc.pp_h4 > 0.1 ? 3 : 2)}
              </p>
              <div className="mt-4 space-y-2 text-[9px]">
                <div className="flex justify-between"><span className="text-[var(--ink-faint)]">PP(H3)</span><span className="font-mono">{pct(selectedColoc.pp_h3, 2)}</span></div>
                <div className="flex justify-between"><span className="text-[var(--ink-faint)]">Sensitivity</span><span className="font-mono">{pct(selectedColoc.pp_h4_sensitivity_min, 2)}–{pct(selectedColoc.pp_h4_sensitivity_max, 2)}</span></div>
                <div className="flex justify-between"><span className="text-[var(--ink-faint)]">Top variant</span><span className="font-mono">{selectedColoc.top_shared_variant}</span></div>
                <div className="flex justify-between"><span className="text-[var(--ink-faint)]">Shared 95% set</span><span className="font-mono">{selectedColoc.shared_credible_set_size_95}</span></div>
              </div>
            </div>
          </div>
          <div className="border-t border-[var(--rule)] bg-[var(--coral-wash)] px-4 py-3">
            <p className="flex gap-2 text-[9px] leading-4 text-[var(--ink-muted)]">
              <TriangleAlert className="mt-0.5 size-3 shrink-0 text-[var(--coral-deep)]" />
              Liver eQTL evidence is near-certain for all three genes, while SORT1 plasma pQTL evidence is weak (PP(H4) 1.64%). Colocalization is evidence, not a naming contest.
            </p>
          </div>
        </div>
      </div>

      <div className="mt-4 grid gap-4 xl:grid-cols-[1fr_320px]">
        <div className="overflow-hidden border border-[var(--rule)] bg-[var(--paper)]">
          <div className="flex items-center justify-between border-b border-[var(--rule)] px-4 py-3">
            <div>
              <p className="text-[10px] font-semibold uppercase tracking-[0.13em]">Gene-level evidence matrix</p>
              <p className="mt-1 text-[9px] text-[var(--ink-faint)]">BayesV2G ranks CELSR2 above SORT1 and PSRC1 in this candidate set</p>
            </div>
            <Braces className="size-4 text-[var(--ink-faint)]" />
          </div>
          <div className="overflow-x-auto p-4">
            <table className="w-full min-w-[690px] border-collapse text-left">
              <thead className="border-b border-[var(--rule)] text-[8px] uppercase tracking-[0.1em] text-[var(--ink-faint)]">
                <tr><th className="pb-2 pr-3 font-medium">Gene</th><th className="pb-2 pr-2 font-medium">Posterior</th><th className="pb-2 px-1 font-medium">eQTL H4</th><th className="pb-2 px-1 font-medium">pQTL H4</th><th className="pb-2 px-1 font-medium">sQTL H4</th><th className="pb-2 px-1 font-medium">VEP</th><th className="pb-2 px-1 font-medium">E2G</th><th className="pb-2 pl-3 text-right font-medium">OT L2G</th></tr>
              </thead>
              <tbody>
                {[...genes].sort((a, b) => a.rank - b.rank).map((gene) => (
                  <tr key={gene.gene_symbol} className="border-b border-[var(--rule-soft)] last:border-b-0">
                    <td className="py-3 pr-3"><span className="text-[10px] font-semibold">{gene.gene_symbol}</span><span className="ml-2 font-mono text-[8px] text-[var(--ink-faint)]">#{gene.rank}</span></td>
                    <td className="py-3 pr-2 font-mono text-[10px] font-semibold">{pct(gene.posterior_mean, 2)}</td>
                    <td className="px-1 py-3"><EvidenceCell value={gene.eQtlColocH4Maximum} /></td>
                    <td className="px-1 py-3"><EvidenceCell value={gene.pQtlColocH4Maximum} /></td>
                    <td className="px-1 py-3"><EvidenceCell value={gene.sQtlColocH4Maximum} /></td>
                    <td className="px-1 py-3"><EvidenceCell value={gene.vepMean} /></td>
                    <td className="px-1 py-3"><EvidenceCell value={gene.e2gMean} /></td>
                    <td className="py-3 pl-3 text-right font-mono text-[9px]">{gene.l2g_score_external.toFixed(3)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>

        <div className="border border-[var(--rule)] bg-[var(--paper)]">
          <div className="border-b border-[var(--rule)] px-4 py-3">
            <p className="text-[10px] font-semibold uppercase tracking-[0.13em]">Sensitivity stability</p>
            <p className="mt-1 text-[9px] text-[var(--ink-faint)]">Ten refits and evidence ablations</p>
          </div>
          <div className="px-4 py-4">
            <div className="flex items-start gap-3">
              <CheckCircle2 className="mt-0.5 size-4 shrink-0 text-[var(--teal-deep)]" />
              <div>
                <p className="text-sm font-semibold">{allStable ? "Top gene unchanged" : "Top gene changed"}</p>
                <p className="mt-1 text-[10px] leading-4 text-[var(--ink-muted)]">CELSR2 remains rank 1 across all {sensitivity.length} SORT1-locus scenarios.</p>
              </div>
            </div>
            <div className="mt-5 border-y border-[var(--rule)] py-3">
              <p className="text-[8px] uppercase tracking-[0.1em] text-[var(--ink-faint)]">Largest probability movement</p>
              <p className="mt-1 font-serif text-2xl">{pct(maxSensitivityShift, 1)}</p>
            </div>
            <p className="mt-4 text-[9px] leading-4 text-[var(--ink-faint)]">
              Stability does not prove the biological mechanism; it shows that this ranking is not driven by a single fragile modeling choice.
            </p>
          </div>
        </div>
      </div>
    </section>
  );
}
