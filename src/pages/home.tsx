import { useState } from "react";
import { Search, ChevronRight } from "lucide-react";
import Loading from "./loading";
import ResultPage from "./ResultPage";

export default function Home() {
  const [prompt, setPrompt] = useState("");
  const [isGenerating, setIsGenerating] = useState(false);
  const [status, setStatus] = useState<string | null>(null);
  const [results, setResults] = useState<any[] | null>(null);

  const handleGenerate = async () => {
    if (!prompt.trim()) {
      setStatus("Please describe your biological intent.");
      return;
    }

    setIsGenerating(true);
    setStatus("Loading...");
    setResults(null);

    try {
      const res = await fetch("http://4.240.107.18/api/drug_discovery", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          text: `(query(${prompt}))`,
        }),
      });

      const json = await res.json();

      if (!Array.isArray(json.results)) {
        throw new Error("Unexpected response format");
      }

      setResults(json.results);
      setStatus("Generation complete.");
    } catch (e) {
      setStatus("ERROR: " + (e as Error).message);
    } finally {
      setIsGenerating(false);
    }
  };

  return (
    <div className="app">
      <style>{`
        * {
          box-sizing: border-box;
          margin: 0;
          padding: 0;
        }

        body {
          font-family: Inter, system-ui, sans-serif;
          background: radial-gradient(circle at top, #0f172a, #020617);
          color: #fff;
        }

        .app {
          min-height: 100vh;
          width: 100vw;
          display: flex;
          align-items: center;
          justify-content: center;
        }

        .hero {
          width: 100%;
          min-height: 100vh;
          display: flex;
          align-items: center;
          justify-content: center;
          padding: 48px;
        }

        .container {
          max-width: 1100px;
          width: 100%;
          text-align: center;
        }

        .hero-title {
          font-size: clamp(40px, 5vw, 64px);
          font-weight: 800;
          line-height: 1.1;
          margin-bottom: 24px;
        }

        .gradient-text {
          background: linear-gradient(135deg, #38bdf8, #22c55e, #a855f7);
          -webkit-background-clip: text;
          -webkit-text-fill-color: transparent;
        }

        .hero-subtitle {
          font-size: 18px;
          color: #cbd5e1;
          max-width: 700px;
          margin: 0 auto 56px;
        }

        .input-section {
          display: flex;
          justify-content: center;
        }

        .input-container {
          background: rgba(15, 23, 42, 0.85);
          border: 1px solid rgba(56, 189, 248, 0.25);
          border-radius: 22px;
          padding: 6px;
          width: 100%;
          max-width: 850px;
          box-shadow: 0 30px 80px rgba(0, 0, 0, 0.45);
        }

        .input-wrapper {
          display: flex;
          align-items: center;
          gap: 16px;
          padding: 18px 22px;
        }

        .input-icon {
          color: #38bdf8;
          flex-shrink: 0;
        }

        .prompt-input {
          flex: 1;
          background: transparent;
          border: none;
          color: white;
          font-size: 16px;
          outline: none;
          resize: none;
          line-height: 1.5;
        }

        .prompt-input::placeholder {
          color: #94a3b8;
        }

        .generate-button {
          display: flex;
          align-items: center;
          gap: 8px;
          background: linear-gradient(135deg, #38bdf8, #6366f1);
          border: none;
          color: white;
          font-size: 16px;
          font-weight: 600;
          padding: 14px 28px;
          border-radius: 16px;
          cursor: pointer;
          transition: all 0.25s ease;
        }

        .generate-button:hover:not(:disabled) {
          transform: translateY(-2px);
          box-shadow: 0 14px 40px rgba(56, 189, 248, 0.5);
        }

        .generate-button:disabled {
          opacity: 0.6;
          cursor: not-allowed;
        }

        .status {
          margin-top: 20px;
          font-size: 14px;
          color: #38bdf8;
        }
      `}</style>

      {isGenerating ? (
        <Loading />
      ) : (
        <div className="hero">
          <div className="container">
            <h1 className="hero-title">
              Design <span className="gradient-text">Molecules</span> with
              <br />
              Artificial Intelligence
            </h1>

            <p className="hero-subtitle">
              Describe symptoms, diseases, or biological targets — our AI generates
              and analyzes candidate molecules with unprecedented speed and accuracy.
            </p>

            <div className="input-section">
              <div className="input-container">
                <div className="input-wrapper">
                  <div className="input-icon">
                    <Search size={24} />
                  </div>

                  <textarea
                    rows={3}
                    className="prompt-input"
                    placeholder="Describe a biological target, disease symptoms, or desired molecular properties…"
                    value={prompt}
                    onChange={(e) => setPrompt(e.target.value)}
                    disabled={isGenerating}
                  />

                  <button
                    className="generate-button"
                    onClick={handleGenerate}
                    disabled={isGenerating}
                  >
                    Generate Molecules
                    <ChevronRight size={20} />
                  </button>
                </div>
              </div>
            </div>

            {status && <div className="status">{status}</div>}

            {results && <ResultPage results={results} />}
          </div>
        </div>
      )}
    </div>
  );
}