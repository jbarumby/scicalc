// SciCalc Pro — execução de SymPy fora da interface.
import { loadPyodide } from "https://cdn.jsdelivr.net/pyodide/v0.27.7/full/pyodide.mjs";

let runtime;
async function initialize() {
  const pyodide = await loadPyodide({
    indexURL: "https://cdn.jsdelivr.net/pyodide/v0.27.7/full/"
  });
  await pyodide.loadPackage("sympy");
  const source = await fetch(new URL("./scicalc-pro-core.py", import.meta.url));
  if (!source.ok) throw new Error("O motor de cálculo não está disponível.");
  await pyodide.runPythonAsync(await source.text());
  return pyodide;
}

runtime = initialize().then(engine => {
  self.postMessage({ kind: "ready" });
  return engine;
}).catch(error => {
  self.postMessage({ kind: "failed", error: String(error?.message || error) });
  throw error;
});

self.onmessage = async event => {
  const { id, request } = event.data || {};
  if (!Number.isInteger(id) || !request) return;
  try {
    const py = await runtime;
    py.globals.set("scicalc_pro_request_json", JSON.stringify(request));
    const raw = await py.runPythonAsync("scicalc_pro_solve(scicalc_pro_request_json)");
    self.postMessage({ kind: "result", id, payload: JSON.parse(raw) });
  } catch (error) {
    self.postMessage({ kind: "result", id, payload: {
      ok: false, error: String(error?.message || error).slice(0, 450)
    }});
  }
};
