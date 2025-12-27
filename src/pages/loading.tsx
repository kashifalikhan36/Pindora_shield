import { useEffect, useRef, useState } from "react";

type Molecule = {
  x: number;
  y: number;
  vx: number;
  vy: number;
  r: number;
  phase: number;
};

export default function Loading() {
  const canvasRef = useRef<HTMLCanvasElement | null>(null);
  const moleculesRef = useRef<Molecule[]>([]);
  const [status, setStatus] = useState("Initializing cellular processes…");

  const MOLECULE_COUNT = 24;
  const INTERACTION_DISTANCE = 90;
  const BASE_SPEED = 2.2;

  /* ===============================
     Init Molecules
  ================================ */
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    moleculesRef.current = Array.from({ length: MOLECULE_COUNT }, () => ({
      x: Math.random() * canvas.width,
      y: Math.random() * canvas.height,
      vx: (Math.random() - 0.5) * BASE_SPEED,
      vy: (Math.random() - 0.5) * BASE_SPEED,
      r: 2.5 + Math.random() * 3,
      phase: Math.random() * Math.PI * 2,
    }));

    let animationId: number;

    const drawMolecule = (m: Molecule) => {
      const pulse = Math.sin(Date.now() * 0.003 + m.phase) * 0.6;

      ctx.beginPath();
      ctx.arc(m.x, m.y, m.r + pulse, 0, Math.PI * 2);
      ctx.fillStyle = "#79e6ff";
      ctx.shadowBlur = 12;
      ctx.shadowColor = "#79e6ff";
      ctx.fill();
    };

    const drawInteraction = (a: Molecule, b: Molecule, d: number) => {
      ctx.strokeStyle = `rgba(121, 230, 255, ${1 - d / INTERACTION_DISTANCE})`;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(a.x, a.y);
      ctx.lineTo(b.x, b.y);
      ctx.stroke();
    };

    const animate = () => {
      ctx.clearRect(0, 0, canvas.width, canvas.height);

      moleculesRef.current.forEach((m, i) => {
        m.vx += (Math.random() - 0.5) * 0.05;
        m.vy += (Math.random() - 0.5) * 0.05;

        m.x += m.vx;
        m.y += m.vy;

        if (m.x < 0 || m.x > canvas.width) m.vx *= -1;
        if (m.y < 0 || m.y > canvas.height) m.vy *= -1;

        for (let j = i + 1; j < moleculesRef.current.length; j++) {
          const o = moleculesRef.current[j];
          const dx = m.x - o.x;
          const dy = m.y - o.y;
          const dist = Math.sqrt(dx * dx + dy * dy);

          if (dist < INTERACTION_DISTANCE) {
            drawInteraction(m, o, dist);
          }
        }

        drawMolecule(m);
      });

      animationId = requestAnimationFrame(animate);
    };

    animate();

    return () => cancelAnimationFrame(animationId);
  }, []);

  /* ===============================
     Status Updates (API Simulation)
  ================================ */
  useEffect(() => {
    const fetchStatus = async () => {
      try {
        const id = Math.floor(Math.random() * 10) + 1;
        const res = await fetch(
          "https://jsonplaceholder.typicode.com/posts/" + id
        );
        const data = await res.json();
        setStatus("Cellular update: " + data.title);
      } catch {
        setStatus("Monitoring biochemical activity…");
      }
    };

    fetchStatus();
    const interval = setInterval(fetchStatus, 4500);
    return () => clearInterval(interval);
  }, []);

  /* ===============================
     UI
  ================================ */
  return (
    <div className="loader-root">
      <style>{`
        .loader-root {
          height: 100vh;
          width: 100vw;
          background: radial-gradient(circle at center, #07131a, #020509);
          display: flex;
          align-items: center;
          justify-content: center;
          font-family: system-ui, sans-serif;
          color: #d9faff;
        }

        .loader {
          text-align: center;
        }

        canvas {
          width: 320px;
          height: 320px;
        }

        .status {
          margin-top: 14px;
          font-size: 13px;
          letter-spacing: 0.3px;
          opacity: 0.85;
        }
      `}</style>

      <div className="loader">
        <canvas ref={canvasRef} width={320} height={320} />
        <div className="status">{status}</div>
      </div>
    </div>
  );
}