# WaveSpecZZ Project — Changelog

## 2025-11-18
- `WaveSpecZZ_1.0.3-pla-kalman-fast-gpuopt-nodetrend.mq5`
  - Exibição exclusiva via `InpViewMode`: ou só ondas (FFT) ou só feed; buffers contrários são preenchidos com `EMPTY_VALUE`.
  - Removido o rótulo/feed label duplicado e qualquer traço do feed quando ondas estão visíveis.
  - Eliminado código de fallback CPU FFT; processamento permanece somente GPU (barra é ignorada se GPU indisponível).
  - `FeedTrace` redimensionado e preenchido apenas no modo feed; propriedades do plot ajustadas (único plot 9).

## 2025-11-17
- `WaveSpecZZ_1.0.3-pla-kalman-fast-gpuopt-nodetrend.mq5`
  - Adicionado combo `InpFeedData` (PLA / ZigZag / Close) como fonte única de `feed_data`.
  - Removidos todos os fallbacks de feed: se PLA ou ZigZag falham, a barra é ignorada e loga erro explícito; modo inválido também aborta a barra.
  - Removidos fallbacks de tempo para ETA: se `PeriodSeconds` ou `delta` entre barras não forem válidos, ETA é desativada na barra e o erro é logado (sem usar 60s).
  - Mantidos logs claros indicando ausência de fallback.
