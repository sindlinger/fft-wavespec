# WaveSpecZZ Project — Changelog

## 2025-11-17
- `WaveSpecZZ_1.0.3-pla-kalman-fast-gpuopt-nodetrend.mq5`
  - Adicionado combo `InpFeedData` (PLA / ZigZag / Close) como fonte única de `feed_data`.
  - Removidos todos os fallbacks de feed: se PLA ou ZigZag falham, a barra é ignorada e loga erro explícito; modo inválido também aborta a barra.
  - Removidos fallbacks de tempo para ETA: se `PeriodSeconds` ou `delta` entre barras não forem válidos, ETA é desativada na barra e o erro é logado (sem usar 60s).
  - Mantidos logs claros indicando ausência de fallback.

