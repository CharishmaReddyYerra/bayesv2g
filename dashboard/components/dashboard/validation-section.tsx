"use client";

import { useState } from "react";
import { ArrowUpRight, BarChart3, ShieldCheck } from "lucide-react";

type Metric = "top1_accuracy" | "recall_at_3" | "mean_reciprocal_rank";

type ValidationRow = {
  model: string;
  n_loci: number;
  top1_accuracy: number;
  recall_at_3: number;
  mean_reciprocal_rank: number;
  multiclass_log_loss: number;
  brier_score: number;
  mean_entropy: number;
};

type CalibrationRow = {
  bin: string;
  mean_confidence: number;
  observed_accuracy: number;
  n: number;
};

type Comparison = {
  top1_difference: number;
  ci_lower_95: number;
  ci_upper_95: number;
  improved_loci: number;
  worsened_loci: number;
  tied_loci: number;
  bootstrap_iterations: number;
};

type Coefficient = {
  feature: string;
  posterior_mean: number;
  posterior_lower: number;
  posterior_upper: number;
  odds_ratio: number;
};

const metricLabels: Record<Metric, string> = {
  top1_accuracy: "Top-1 accuracy",
  recall_at_3: "Recall @ 3",
  mean_reciprocal_rank: "Mean reciprocal rank",
};

const modelLabels: Record<string, string> = {
  BayesV2G: "BayesV2G",
  distance_only: "Distance only",
  uniform: "Uniform prior",
  eqtl_coloc_only: "eQTL coloc only",
  vep_only: "VEP only",
};

const featureLabels: Record<string, string> = {
  distanceTssMean: "TSS proximity",
  eQtlColocH4Maximum: "eQTL coloc H4",
  pQtlColocH4Maximum: "pQTL coloc H4",
  sQtlColocH4Maximum: "sQTL coloc H4",
  vepMean: "VEP consequence",
  e2gMean: "Enhancer-to-gene",
};

function pct(value: number, digits = 1) {
  return `${(value * 100).toFixed(digits)}%`;
}

export function ValidationSection({
  validation,
  calibration,
  comparison,
  coefficients,
}: {
  validation: ValidationRow[];
  calibration: CalibrationRow[];
  comparison: Comparison;
  coefficients: Coefficient[];
}) {
  const [metric, setMetric] = useState<Metric>("top1_accuracy");
  const ranked = [...validation].sort((a, b) => b[metric] - a[metric]);
  const bayes = validation.find((row) => row.model === "BayesV2G")!;
  const maxCoefficient = Math.max(...coefficients.map((row) => row.posterior_upper));
  const x = (value: number) => 28 + value * 166;
  const y = (value: number) => 192 - value * 166;
  const calibrationPath = calibration
    .map((row, index) => `${index === 0 ? "M" : "L"}${x(row.mean_confidence)},${y(row.observed_accuracy)}`)
    .join(" ");

  const ciDomain = { min: -0.03, max: 0.13 };
  const ciX = (value: number) => ((value - ciDomain.min) / (ciDomain.max - ciDomain.min)) * 100;

  return (
    <section id="validation" className="scroll-mt-20 border-t border-[var(--rule-strong)] pt-10">
      <div className="mb-6 grid gap-4 lg:grid-cols-[1fr_0.72fr] lg:items-end">
        <div>
          <p className="mb-2 flex items-center gap-2 text-[10px] font-semibold uppercase tracking-[0.16em] text-[var(--cobalt)]">
            <ShieldCheck className="size-3.5" /> Held-out evidence
          </p>
          <h2 className="font-serif text-3xl tracking-[-0.03em]">Validation without the victory lap.</h2>
        </div>
        <p className="text-xs leading-5 text-[var(--ink-muted)]">
          Performance is reported on 1,048 held-out Open Targets loci. The paired result is bootstrapped at the locus level; calibration uses the winning gene probability.
        </p>
      </div>

      <div className="grid gap-4 xl:grid-cols-[1.05fr_0.9fr_0.82fr]">
        <div className="border border-[var(--rule)] bg-[var(--paper)]">
          <div className="flex flex-wrap items-center justify-between gap-3 border-b border-[var(--rule)] px-4 py-3">
            <div>
              <p className="text-[10px] font-semibold uppercase tracking-[0.13em]">Benchmark models</p>
              <p className="mt-1 text-[9px] text-[var(--ink-faint)]">Same 1,048-locus test partition</p>
            </div>
            <div className="inline-flex border border-[var(--rule)] p-0.5">
              {(Object.keys(metricLabels) as Metric[]).map((key) => (
                <button
                  key={key}
                  type="button"
                  onClick={() => setMetric(key)}
                  className={`px-2 py-1.5 text-[9px] transition-colors ${
                    metric === key ? "bg-[var(--ink)] text-white" : "text-[var(--ink-muted)] hover:bg-[var(--canvas)]"
                  }`}
                >
                  {key === "top1_accuracy" ? "Top-1" : key === "recall_at_3" ? "R@3" : "MRR"}
                </button>
              ))}
            </div>
          </div>
          <div className="space-y-4 p-4">
            {ranked.map((row, index) => (
              <div key={row.model} className="grid grid-cols-[105px_1fr_48px] items-center gap-3">
                <span className={`text-[10px] ${row.model === "BayesV2G" ? "font-semibold text-[var(--ink)]" : "text-[var(--ink-muted)]"}`}>
                  <span className="mr-2 font-mono text-[8px] text-[var(--ink-faint)]">{String(index + 1).padStart(2, "0")}</span>
                  {modelLabels[row.model] ?? row.model}
                </span>
                <span className="relative h-2 bg-[var(--rule-soft)]">
                  <span
                    className={`absolute inset-y-0 left-0 ${row.model === "BayesV2G" ? "bg-[var(--cobalt)]" : "bg-[var(--slate)]"}`}
                    style={{ width: `${row[metric] * 100}%` }}
                  />
                </span>
                <span className="text-right font-mono text-[10px] font-semibold tabular-nums">{pct(row[metric])}</span>
              </div>
            ))}
          </div>
          <div className="grid grid-cols-3 gap-px border-t border-[var(--rule)] bg-[var(--rule)]">
            <div className="bg-[var(--canvas)] px-3 py-3">
              <p className="text-[8px] uppercase tracking-[0.1em] text-[var(--ink-faint)]">BayesV2G top-1</p>
              <p className="mt-1 font-mono text-base font-semibold">{pct(bayes.top1_accuracy, 2)}</p>
            </div>
            <div className="bg-[var(--canvas)] px-3 py-3">
              <p className="text-[8px] uppercase tracking-[0.1em] text-[var(--ink-faint)]">Recall @ 3</p>
              <p className="mt-1 font-mono text-base font-semibold">{pct(bayes.recall_at_3, 2)}</p>
            </div>
            <div className="bg-[var(--canvas)] px-3 py-3">
              <p className="text-[8px] uppercase tracking-[0.1em] text-[var(--ink-faint)]">MRR</p>
              <p className="mt-1 font-mono text-base font-semibold">{pct(bayes.mean_reciprocal_rank, 2)}</p>
            </div>
          </div>
        </div>

        <div className="border border-[var(--rule)] bg-[var(--paper)]">
          <div className="border-b border-[var(--rule)] px-4 py-3">
            <p className="text-[10px] font-semibold uppercase tracking-[0.13em]">Reliability</p>
            <p className="mt-1 text-[9px] text-[var(--ink-faint)]">Confidence versus observed top-gene accuracy</p>
          </div>
          <div className="p-3">
            <svg viewBox="0 0 222 215" className="mx-auto block w-full max-w-[310px]" role="img" aria-label="BayesV2G calibration curve">
              <line x1="28" y1="192" x2="194" y2="26" stroke="var(--rule-strong)" strokeDasharray="4 4" />
              {[0, 0.5, 1].map((value) => (
                <g key={value}>
                  <line x1={x(value)} y1="22" x2={x(value)} y2="192" stroke="var(--rule-soft)" />
                  <line x1="28" y1={y(value)} x2="198" y2={y(value)} stroke="var(--rule-soft)" />
                  <text x={x(value)} y="208" textAnchor="middle" fill="var(--ink-faint)" fontSize="8" fontFamily="monospace">{value.toFixed(1)}</text>
                  <text x="20" y={y(value) + 3} textAnchor="end" fill="var(--ink-faint)" fontSize="8" fontFamily="monospace">{value.toFixed(1)}</text>
                </g>
              ))}
              <path d={calibrationPath} fill="none" stroke="var(--cobalt)" strokeWidth="2" />
              {calibration.map((row) => (
                <circle
                  key={row.bin}
                  cx={x(row.mean_confidence)}
                  cy={y(row.observed_accuracy)}
                  r={Math.min(3 + Math.sqrt(row.n) / 5, 8)}
                  fill="var(--paper)"
                  stroke="var(--cobalt)"
                  strokeWidth="1.5"
                >
                  <title>{`${row.bin}: confidence ${pct(row.mean_confidence)}, observed ${pct(row.observed_accuracy)}, n=${row.n}`}</title>
                </circle>
              ))}
              <text x="111" y="214" textAnchor="middle" fill="var(--ink-muted)" fontSize="8">mean confidence</text>
              <text transform="translate(8 110) rotate(-90)" textAnchor="middle" fill="var(--ink-muted)" fontSize="8">observed accuracy</text>
            </svg>
          </div>
        </div>

        <div className="border border-[var(--rule)] bg-[var(--paper)]">
          <div className="border-b border-[var(--rule)] px-4 py-3">
            <p className="text-[10px] font-semibold uppercase tracking-[0.13em]">Paired top-1 effect</p>
            <p className="mt-1 text-[9px] text-[var(--ink-faint)]">BayesV2G minus distance-only</p>
          </div>
          <div className="px-4 py-5">
            <p className="flex items-center gap-1 font-serif text-4xl tracking-[-0.04em] text-[var(--teal-deep)]">
              +{pct(comparison.top1_difference, 2)} <ArrowUpRight className="size-4" />
            </p>
            <p className="mt-1 text-[10px] text-[var(--ink-muted)]">
              95% CI {pct(comparison.ci_lower_95, 2)} to {pct(comparison.ci_upper_95, 2)}
            </p>
            <div className="relative mt-7 h-10 border-y border-[var(--rule-soft)]">
              <span className="absolute inset-y-0 w-px bg-[var(--rule-strong)]" style={{ left: `${ciX(0)}%` }} />
              <span
                className="absolute top-1/2 h-px -translate-y-1/2 bg-[var(--teal-deep)]"
                style={{ left: `${ciX(comparison.ci_lower_95)}%`, width: `${ciX(comparison.ci_upper_95) - ciX(comparison.ci_lower_95)}%` }}
              />
              <span
                className="absolute top-1/2 size-2.5 -translate-x-1/2 -translate-y-1/2 rotate-45 bg-[var(--teal-deep)]"
                style={{ left: `${ciX(comparison.top1_difference)}%` }}
              />
            </div>
            <div className="mt-2 flex justify-between font-mono text-[8px] text-[var(--ink-faint)]">
              <span>−3 pp</span><span>0</span><span>+13 pp</span>
            </div>
            <div className="mt-6 grid grid-cols-3 gap-px border border-[var(--rule)] bg-[var(--rule)]">
              <div className="bg-[var(--canvas)] px-2 py-2.5 text-center">
                <p className="font-mono text-sm font-semibold text-[var(--teal-deep)]">{comparison.improved_loci}</p>
                <p className="mt-0.5 text-[8px] text-[var(--ink-faint)]">improved</p>
              </div>
              <div className="bg-[var(--canvas)] px-2 py-2.5 text-center">
                <p className="font-mono text-sm font-semibold text-[var(--coral-deep)]">{comparison.worsened_loci}</p>
                <p className="mt-0.5 text-[8px] text-[var(--ink-faint)]">worsened</p>
              </div>
              <div className="bg-[var(--canvas)] px-2 py-2.5 text-center">
                <p className="font-mono text-sm font-semibold">{comparison.tied_loci}</p>
                <p className="mt-0.5 text-[8px] text-[var(--ink-faint)]">tied</p>
              </div>
            </div>
            <p className="mt-3 text-[8px] leading-4 text-[var(--ink-faint)]">
              {comparison.bootstrap_iterations.toLocaleString()} paired locus bootstrap iterations.
            </p>
          </div>
        </div>
      </div>

      <div className="mt-4 border border-[var(--rule)] bg-[var(--paper)]">
        <div className="grid border-b border-[var(--rule)] px-4 py-3 lg:grid-cols-[220px_1fr] lg:items-center">
          <div className="flex items-center gap-2">
            <BarChart3 className="size-3.5 text-[var(--cobalt)]" />
            <p className="text-[10px] font-semibold uppercase tracking-[0.13em]">Posterior coefficients</p>
          </div>
          <p className="mt-1 text-[9px] text-[var(--ink-faint)] lg:mt-0">Predictive log-odds weights; magnitudes are not causal effect sizes.</p>
        </div>
        <div className="grid gap-x-6 gap-y-3 p-4 sm:grid-cols-2 xl:grid-cols-3">
          {coefficients.map((row) => (
            <div key={row.feature} className="grid grid-cols-[105px_1fr_40px] items-center gap-2">
              <span className="text-[9px] text-[var(--ink-muted)]">{featureLabels[row.feature] ?? row.feature}</span>
              <span className="relative h-1.5 bg-[var(--rule-soft)]">
                <span className="absolute inset-y-0 left-0 bg-[var(--cobalt)]" style={{ width: `${(row.posterior_mean / maxCoefficient) * 100}%` }} />
                <span className="absolute top-1/2 h-3 -translate-y-1/2 border-x border-[var(--ink-muted)]" style={{ left: `${(row.posterior_lower / maxCoefficient) * 100}%`, width: `${((row.posterior_upper - row.posterior_lower) / maxCoefficient) * 100}%` }} />
              </span>
              <span className="text-right font-mono text-[9px] font-semibold">{row.posterior_mean.toFixed(2)}</span>
            </div>
          ))}
        </div>
      </div>
    </section>
  );
}
