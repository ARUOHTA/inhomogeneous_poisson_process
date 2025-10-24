# Agent Notes

- When inspecting files/notebooks, avoid using `python - <<'PY' ...` unless necessary.
- Prefer fast tools (`rg`, `sed`, `head`, `tail`) for quick lookups to minimize overhead.
- Before editing, confirm contexts with `sed`/`rg` to avoid heavy Python invocations unless required.
