"use client";

import { useMemo, useState } from "react";

type Variant = {
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
};

function pct(value: number) {
  if (value >= 0.01) return `${(value * 100).toFixed(2)}%`;
  if (value >= 0.0001) return `${(value * 100).toFixed(3)}%`;
  return value.toExponential(1);
}

export function CredibleVariantTrack({
  variants,
  leadVariant,
}: {
  variants: Variant[];
  leadVariant: string;
}) {
  const ordered = useMemo(
    () => [...variants].sort((a, b) => a.position - b.position),
    [variants],
  );
  const ranked = useMemo(
    () => [...variants].sort((a, b) => b.posterior_inclusion_probability - a.posterior_inclusion_probability),
    [variants],
  );
  const [selectedId, setSelectedId] = useState(leadVariant);
  const selected = ranked.find((variant) => variant.rsid === selectedId) ?? ranked[0];
  const min = ordered[0]?.position ?? 0;
  const max = ordered.at(-1)?.position ?? min + 1;
  const x = (position: number) => 2 + ((position - min) / Math.max(max - min, 1)) * 96;

  return (
    <div className="mt-5 border-t border-[var(--rule)] pt-4">
      <div className="flex flex-wrap items-end justify-between gap-3">
        <div>
          <p className="text-[10px] font-semibold uppercase tracking-[0.14em] text-[var(--ink-muted)]">
            Fine-mapped variant architecture
          </p>
          <p className="mt-1 text-[10px] text-[var(--ink-faint)]">
            Posterior inclusion probability across {variants.length} displayed credible-set variants
          </p>
        </div>
        {selected && (
          <p className="text-right font-mono text-[9px] leading-4 text-[var(--ink-faint)]">
            <strong className="font-semibold text-[var(--ink)]">{selected.rsid}</strong> · {selected.reference_allele}›{selected.alternate_allele}
            <br />PIP {pct(selected.posterior_inclusion_probability)} · {selected.most_severe_consequence.replaceAll("_", " ")}
          </p>
        )}
      </div>

      <div className="relative mt-4 h-[82px] border-b border-[var(--rule-strong)]" role="img" aria-label="Fine-mapped variant posterior inclusion probabilities">
        <div className="absolute inset-x-0 top-[20%] border-t border-dashed border-[var(--rule-soft)]" />
        <div className="absolute inset-x-0 top-[50%] border-t border-dashed border-[var(--rule-soft)]" />
        {ordered.map((variant) => {
          const active = variant.rsid === selected.rsid;
          return (
            <button
              key={`${variant.rsid}-${variant.position}`}
              type="button"
              onClick={() => setSelectedId(variant.rsid)}
              aria-label={`${variant.rsid}, PIP ${pct(variant.posterior_inclusion_probability)}`}
              className={`absolute bottom-0 w-[3px] min-w-[3px] -translate-x-1/2 transition-colors ${
                active ? "z-10 bg-[var(--coral)]" : variant.is_95_credible_set ? "bg-[var(--cobalt)]" : "bg-[var(--slate)]"
              }`}
              style={{
                left: `${x(variant.position)}%`,
                height: `${Math.max(variant.posterior_inclusion_probability * 100, 3)}%`,
              }}
            />
          );
        })}
      </div>
      <div className="mt-2 flex justify-between font-mono text-[8px] text-[var(--ink-faint)]">
        <span>{min.toLocaleString()}</span>
        <span>{max.toLocaleString()}</span>
      </div>
    </div>
  );
}
