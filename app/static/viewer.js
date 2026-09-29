"use strict";
// Original, dependency-free projection of supplied molecular coordinates.
// This renderer does not generate geometries or infer chemical properties.
class MoleculeViewer {
  static colors = { H: "#d8dfd0", C: "#64805b", N: "#628bbb", O: "#cb7565", F: "#9abe69", P: "#ba9860", S: "#d3bd66", Cl: "#76a26d", Br: "#ad7961", I: "#9b7dab" };
  constructor(canvas, legend) {
    this.canvas = canvas; this.legend = legend; this.ctx = canvas.getContext("2d");
    this.atoms = []; this.bonds = []; this.points = []; this.yaw = .45; this.pitch = -.3; this.zoom = 1; this.frame = null; this.pointers = new Map(); this.gestureDistance = null;
    this.resizeObserver = new ResizeObserver(() => this.render()); this.resizeObserver.observe(canvas);
    canvas.addEventListener("pointerdown", (event) => { canvas.focus({ preventScroll: true }); canvas.setPointerCapture(event.pointerId); this.pointers.set(event.pointerId, [event.clientX, event.clientY]); this.gestureDistance = this.pointerDistance(); });
    canvas.addEventListener("pointermove", (event) => {
      const previous = this.pointers.get(event.pointerId); if (!previous) return;
      this.pointers.set(event.pointerId, [event.clientX, event.clientY]);
      if (this.pointers.size > 1) { const distance = this.pointerDistance(); if (this.gestureDistance > 0) this.zoom = Math.min(4, Math.max(.45, this.zoom * distance / this.gestureDistance)); this.gestureDistance = distance; }
      else { this.yaw += (event.clientX - previous[0]) * .012; this.pitch += (event.clientY - previous[1]) * .012; }
      this.render();
    });
    const release = (event) => { this.pointers.delete(event.pointerId); this.gestureDistance = this.pointerDistance(); };
    canvas.addEventListener("pointerup", release); canvas.addEventListener("pointercancel", release); canvas.addEventListener("lostpointercapture", release);
    canvas.addEventListener("wheel", (event) => { event.preventDefault(); this.scale(Math.exp(-event.deltaY * .0015)); }, { passive: false });
    canvas.addEventListener("keydown", (event) => {
      const actions = { ArrowLeft: () => { this.yaw -= .15; }, ArrowRight: () => { this.yaw += .15; }, ArrowUp: () => { this.pitch -= .15; }, ArrowDown: () => { this.pitch += .15; }, "+": () => this.scale(1.12), "=": () => this.scale(1.12), "-": () => this.scale(.89), r: () => this.reset(), R: () => this.reset() };
      if (actions[event.key]) { event.preventDefault(); actions[event.key](); this.render(); }
    });
    document.addEventListener("themechange", () => this.render());
  }
  pointerDistance() { if (this.pointers.size < 2) return null; const [a, b] = [...this.pointers.values()]; return Math.hypot(a[0] - b[0], a[1] - b[1]); }
  setStructure(atoms, bonds, coordinates) {
    if (!Array.isArray(atoms) || !Array.isArray(coordinates) || atoms.length !== coordinates.length || coordinates.some((p) => !Array.isArray(p) || p.length !== 3 || !p.every(Number.isFinite))) throw new Error("三维坐标不完整，无法显示。");
    this.atoms = atoms; this.bonds = Array.isArray(bonds) ? bonds : [];
    const center = [0, 1, 2].map((axis) => coordinates.reduce((sum, p) => sum + p[axis], 0) / Math.max(1, coordinates.length));
    this.points = coordinates.map((p) => p.map((v, axis) => v - center[axis]));
    this.radius = Math.max(1, ...this.points.map((p) => Math.hypot(...p)));
    this.legend.replaceChildren();
    [...new Set(atoms.map((a) => a.element))].forEach((element) => { const label = document.createElement("span"), dot = document.createElement("i"); dot.style.background = MoleculeViewer.colors[element] || "#af93a4"; label.append(dot, document.createTextNode(element)); this.legend.append(label); });
    this.render();
  }
  reset() { this.yaw = .45; this.pitch = -.3; this.zoom = 1; this.render(); }
  scale(factor) { this.zoom = Math.max(.45, Math.min(4, this.zoom * factor)); this.render(); }
  render() { if (this.frame !== null) return; this.frame = requestAnimationFrame(() => { this.frame = null; this.draw(); }); }
  draw() {
    const rect = this.canvas.getBoundingClientRect(); if (!rect.width || !rect.height || !this.ctx) return;
    const ratio = Math.min(window.devicePixelRatio || 1, 2), width = rect.width, height = rect.height;
    if (this.canvas.width !== Math.round(width * ratio) || this.canvas.height !== Math.round(height * ratio)) { this.canvas.width = Math.round(width * ratio); this.canvas.height = Math.round(height * ratio); }
    const ctx = this.ctx; ctx.setTransform(ratio, 0, 0, ratio, 0, 0); ctx.clearRect(0, 0, width, height); if (!this.points.length) return;
    const fit = Math.min(width - 65, height - 60) / (this.radius * 2.35) * this.zoom, cy = Math.cos(this.yaw), sy = Math.sin(this.yaw), cx = Math.cos(this.pitch), sx = Math.sin(this.pitch);
    const projected = this.points.map(([x, y, z], i) => { const rx = x * cy + z * sy, rz = -x * sy + z * cy, ry = y * cx - rz * sx, depth = y * sx + rz * cx, perspective = 1 / (1 - depth / (this.radius * 7)); return { i, x: width / 2 + rx * fit * perspective, y: height / 2 - ry * fit * perspective, z: depth, r: Math.max(2.5, Math.min(18, (this.atoms[i].element === "H" ? .16 : .28) * fit * perspective)) }; });
    const shapes = projected.map((p) => ({ type: "atom", z: p.z, p }));
    this.bonds.forEach((bond) => { const a = projected[bond.begin], b = projected[bond.end]; if (a && b) shapes.push({ type: "bond", z: (a.z + b.z) / 2 - .12, a, b, order: bond.order }); });
    shapes.sort((a, b) => a.z - b.z);
    const dark = document.documentElement.dataset.theme === "dark";
    shapes.forEach((shape) => {
      if (shape.type === "bond") {
        const { a, b } = shape, length = Math.hypot(b.x - a.x, b.y - a.y) || 1, nx = -(b.y - a.y) / length, ny = (b.x - a.x) / length;
        const n = shape.order >= 2.8 ? 3 : shape.order >= 1.8 ? 2 : 1;
        ctx.lineCap = "round"; ctx.lineWidth = Math.max(1.3, Math.min(4.8, fit * .08));
        const gradient = ctx.createLinearGradient(a.x, a.y, b.x, b.y); gradient.addColorStop(0, MoleculeViewer.colors[this.atoms[a.i].element] || "#829877"); gradient.addColorStop(1, MoleculeViewer.colors[this.atoms[b.i].element] || "#829877"); ctx.strokeStyle = gradient;
        for (let line = 0; line < n; line++) { const offset = (line - (n - 1) / 2) * Math.max(2.5, fit * .1); ctx.beginPath(); ctx.moveTo(a.x + nx * offset, a.y + ny * offset); ctx.lineTo(b.x + nx * offset, b.y + ny * offset); ctx.stroke(); }
      } else {
        const { p } = shape, element = this.atoms[p.i].element, color = MoleculeViewer.colors[element] || "#a48ba5";
        const gradient = ctx.createRadialGradient(p.x - p.r * .3, p.y - p.r * .35, p.r * .05, p.x, p.y, p.r); gradient.addColorStop(0, dark ? "#e3e9db" : "#f9fcf6"); gradient.addColorStop(.35, color); gradient.addColorStop(1, color);
        ctx.beginPath(); ctx.arc(p.x, p.y, p.r, 0, Math.PI * 2); ctx.fillStyle = gradient; ctx.fill(); ctx.strokeStyle = dark ? "#1c281dbb" : "#405b3830"; ctx.lineWidth = .7; ctx.stroke();
        if (element !== "C" && element !== "H" && p.r >= 5) { ctx.fillStyle = "#fff"; ctx.font = `600 ${Math.max(8, Math.min(12, p.r))}px Segoe UI`; ctx.textAlign = "center"; ctx.textBaseline = "middle"; ctx.fillText(element, p.x, p.y + .5); }
      }
    });
  }
}
window.MoleculeViewer = MoleculeViewer;
