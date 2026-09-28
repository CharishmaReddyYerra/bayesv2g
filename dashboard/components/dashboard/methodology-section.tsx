import { Archive, Binary, Database, FileCheck2, Scale } from "lucide-react";

const sources = [
  { name: "Open Targets L2G 26.06", detail: "99,048 training · 16,804 test gene–locus pairs" },
  { name: "Open Targets 26.09", detail: "10 LDL loci · 26 surfaced candidates" },
  { name: "GWAS Catalog", detail: "GCST90018961 · 5,772 SORT1-region variants" },
  { name: "eQTL Catalogue", detail: "GTEx liver 8,631 · plasma pQTL 3,892 rows" },
  { name: "Ensembl", detail: "Gene annotation, VEP consequence, and LD services" },
];

const caveats = [
  ["Candidate ascertainment", "Genes are externally surfaced at Open Targets L2G ≥ 0.05. Omitted genes cannot rank."],
  ["Unresolved mass", "Absolute probabilities reserve 20% prior mass and account for credible-set coverage."],
  ["Singleton sets", "APOB and LDLR have one surfaced candidate; conditional probability 1.0 is structural."],
  ["Evidence zeros", "A zero may mean no signal, missing tissue, or low power—not biological exclusion."],
  ["Mechanistic scope", "The model prioritizes a primary target distribution even when loci may have multiple mediators."],
  ["LD and coloc", "Fine-mapping uses out-of-sample LD; the independent coloc audit assumes one causal variant."],
];

export function MethodologySection() {
  return (
    <section id="methods" className="scroll-mt-20 border-t border-[var(--rule-strong)] pt-10">
      <div className="mb-6 grid gap-4 lg:grid-cols-[1fr_0.72fr] lg:items-end">
        <div>
          <p className="mb-2 flex items-center gap-2 text-[10px] font-semibold uppercase tracking-[0.16em] text-[var(--teal-deep)]">
            <FileCheck2 className="size-3.5" /> Provenance and limits
          </p>
          <h2 className="font-serif text-3xl tracking-[-0.03em]">What the number means—and what it does not.</h2>
        </div>
        <p className="text-xs leading-5 text-[var(--ink-muted)]">
          The dashboard reads the checked-in analysis outputs directly. The website adds interaction and interpretation; it does not recompute or alter the posterior.
        </p>
      </div>

      <div className="grid gap-4 xl:grid-cols-[0.8fr_1.2fr]">
        <div className="border border-[var(--rule)] bg-[var(--ink)] text-white">
          <div className="border-b border-white/15 px-5 py-4">
            <p className="flex items-center gap-2 text-[10px] font-semibold uppercase tracking-[0.13em] text-white/65">
              <Binary className="size-3.5" /> Probability construction
            </p>
          </div>
          <div className="px-5 py-6">
            <p className="font-serif text-[22px] leading-[1.45] tracking-[-0.02em]">
              P<sub className="text-xs">absolute</sub> = 0.80 × credible-set coverage × P<sub className="text-xs">conditional</sub>
            </p>
            <div className="my-5 h-px bg-white/15" />
            <p className="text-[10px] leading-5 text-white/65">
              The residual is explicit unresolved probability. This keeps a strong within-candidate ranking from masquerading as certainty that the true effector is present in the surfaced set.
            </p>
          </div>
          <div className="grid grid-cols-3 gap-px bg-white/15">
            <div className="bg-[var(--ink)] px-3 py-3"><p className="font-mono text-base font-semibold">10</p><p className="mt-1 text-[8px] text-white/50">LDL loci</p></div>
            <div className="bg-[var(--ink)] px-3 py-3"><p className="font-mono text-base font-semibold">26</p><p className="mt-1 text-[8px] text-white/50">gene hypotheses</p></div>
            <div className="bg-[var(--ink)] px-3 py-3"><p className="font-mono text-base font-semibold">120</p><p className="mt-1 text-[8px] text-white/50">credible variants</p></div>
          </div>
        </div>

        <div className="border border-[var(--rule)] bg-[var(--paper)]">
          <div className="flex items-center justify-between border-b border-[var(--rule)] px-4 py-3">
            <div>
              <p className="text-[10px] font-semibold uppercase tracking-[0.13em]">Source ledger</p>
              <p className="mt-1 text-[9px] text-[var(--ink-faint)]">Pinned real-data snapshots used by the reproducible run</p>
            </div>
            <Database className="size-4 text-[var(--ink-faint)]" />
          </div>
          <div>
            {sources.map((source, index) => (
              <div key={source.name} className="grid grid-cols-[26px_1fr] gap-3 border-b border-[var(--rule-soft)] px-4 py-3 last:border-b-0 sm:grid-cols-[26px_180px_1fr]">
                <span className="font-mono text-[9px] text-[var(--ink-faint)]">{String(index + 1).padStart(2, "0")}</span>
                <span className="text-[10px] font-semibold">{source.name}</span>
                <span className="col-start-2 text-[9px] text-[var(--ink-muted)] sm:col-start-auto">{source.detail}</span>
              </div>
            ))}
          </div>
          <div className="flex items-center gap-2 border-t border-[var(--rule)] bg-[var(--canvas)] px-4 py-3 text-[9px] text-[var(--ink-muted)]">
            <Archive className="size-3.5" /> Run manifest, exact release IDs, and SHA-256 checksums are included with the analysis outputs.
          </div>
        </div>
      </div>

      <div className="mt-4 border border-[var(--rule)] bg-[var(--paper)]">
        <div className="flex items-center gap-2 border-b border-[var(--rule)] px-4 py-3">
          <Scale className="size-3.5 text-[var(--coral-deep)]" />
          <p className="text-[10px] font-semibold uppercase tracking-[0.13em]">Interpretive guardrails</p>
        </div>
        <div className="grid md:grid-cols-2 xl:grid-cols-3">
          {caveats.map(([title, detail], index) => (
            <div key={title} className={`px-4 py-4 ${index < caveats.length - 1 ? "border-b border-[var(--rule-soft)] xl:border-b-0" : ""} ${index % 3 !== 2 ? "xl:border-r" : ""} ${index % 2 === 0 ? "md:border-r xl:border-r" : "md:border-r-0"}`}>
              <p className="text-[10px] font-semibold">{title}</p>
              <p className="mt-1.5 text-[9px] leading-4 text-[var(--ink-muted)]">{detail}</p>
            </div>
          ))}
        </div>
      </div>

      <footer className="mt-8 flex flex-col justify-between gap-3 border-t border-[var(--rule)] py-5 text-[9px] text-[var(--ink-faint)] sm:flex-row sm:items-center">
        <span>BayesV2G Evidence Workbench · LDL demonstration · reproducible public-data snapshot</span>
        <span className="font-semibold uppercase tracking-[0.11em] text-[var(--coral-deep)]">Research use only · not clinical guidance</span>
      </footer>
    </section>
  );
}
