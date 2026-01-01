import { useEffect, useRef, useState } from "react";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";




type Molecule = {
  smiles: string;
  similarity: number;
  properties: {
    molecular_weight: number;
    logp: number;
    hbd: number;
    hba: number;
    rotatable_bonds: number;
    aromatic_rings: number;
  };
};

type ResultItem = {
  disease_name: string;
  target_symbol: string;
  drug_name: string;
  generated_molecules: Molecule[][];
};

export default function ResultPage({ results }: { results: ResultItem[] }) {
  const topRef = useRef<HTMLDivElement | null>(null);

  const [metrics, setMetrics] = useState<any | null>(null);
  const [activeSmile, setActiveSmile] = useState<string | null>(null);
  const [isMetricsLoading, setIsMetricsLoading] = useState(false);


  const fetchMetrics = async (smile: string) => {
  try {
    setActiveSmile(smile);        // jis SMILES pe click hua
    setMetrics(null);             // purana data hatao
    setIsMetricsLoading(true);

    const res = await fetch(
      "http://4.240.107.18/metrics/metrics_data"
,
      {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ input_smile: smile }),
      }
    );

    const data = await res.json();
    setMetrics(data);
  } catch (err) {
    console.error("Metrics API error:", err);
  } finally {
    setIsMetricsLoading(false);
  }
};



  useEffect(() => {
    if (topRef.current) {
      topRef.current.scrollIntoView({
        behavior: "smooth",
        block: "start",
      });
    }
  }, [results]);

  return (
    <div ref={topRef} className="results-page">
      <style>{`
        .results-page {
          margin-top: 72px;
          display: grid;
          gap: 32px;
          max-width: 1100px;
          width: 100%;
          margin-left: auto;
          margin-right: auto;
        }

        .result-card {
          background: rgba(15, 23, 42, 0.9);
          border: 1px solid rgba(148, 163, 184, 0.15);
          border-radius: 22px;
          padding: 28px;
          text-align: left;
        }

        .result-header h3 {
          font-size: 22px;
          margin-bottom: 6px;
          font-weight: 700;
        }

        .result-meta {
          font-size: 14px;
          color: #94a3b8;
          margin-bottom: 18px;
          line-height: 1.6;
        }

        .props {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
          gap: 12px;
          font-size: 14px;
          color: #e5e7eb;
          background: rgba(2, 6, 23, 0.6);
          border: 1px solid rgba(148, 163, 184, 0.12);
          padding: 16px;
          border-radius: 14px;
        }

        .smiles {
          margin-top: 14px;
          font-family: monospace;
          font-size: 13px;
          color: #38bdf8;
          word-break: break-all;
        }
      `}</style>

      {results.map((item, index) => {
        const mol = item.generated_molecules?.[0]?.[0];

        return (
          <div key={index} className="result-card">
            <div className="result-header">
              <h3>
                {index + 1}. {item.drug_name}
              </h3>
            </div>

            <div className="result-meta">
              Disease: {item.disease_name}
              <br />
              Target: {item.target_symbol}
            </div>

            {mol && (
              <>
                <div className="props">
                  <div>Molecular Weight: {mol.properties.molecular_weight}</div>
                  <div>LogP: {mol.properties.logp}</div>
                  <div>HBD: {mol.properties.hbd}</div>
                  <div>HBA: {mol.properties.hba}</div>
                  <div>Rotatable Bonds: {mol.properties.rotatable_bonds}</div>
                  <div>Aromatic Rings: {mol.properties.aromatic_rings}</div>
                </div>

                <div
                  className="smiles"
                  style={{ cursor: "pointer" }}
                  onClick={() => fetchMetrics(mol.smiles)}
                >
                  SMILES: {mol.smiles}
                </div>
                
                {activeSmile === mol.smiles && (
                  <div className="props" style={{ marginTop: "16px" }}>
                    {isMetricsLoading ? (
                      <div>Analyzing selected molecule…</div>
                    ) : metrics ? (
                      <div className="prose prose-invert max-w-none text-sm">
                        <ReactMarkdown
                          remarkPlugins={[remarkMath]}
                          rehypePlugins={[rehypeKatex]}
                        >
                          {metrics?.report}
                        </ReactMarkdown>
                      </div>


                    ) : null}
                  </div>
                )}

              </>
            )}
          </div>
        );
      })}
    </div>
  );
}
