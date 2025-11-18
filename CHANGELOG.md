# WaveSpecZZ Project — Changelog

## 2025-11-18
- `WaveSpecZZ_1.1.0-gpuopt.mq5`
  - Import do bridge atualizado para `gpu_extract_cycles`; o indicador agora obtém os ciclos diretamente na GPU (FFT ou MUSIC/ESPRIT) conforme `InpGpuMethod/InpGpuArOrder`, sem qualquer fallback CPU.
  - Defaults ajustados para “ciclos perfeitos”: `InpGpuMethod=1` (MUSIC/ESPRIT), `InpGpuArOrder=10`, `InpGpuTopK=2`, janela 4096.
  - Limpamos buffers por barra antes de desenhar para evitar linhas horizontais quando menos de 8 ciclos vêm do core.
  - Mantida exibição exclusiva (feed ou waves); feed ocultado quando ondas ativas.
- `WaveSpecZZ_1.0.3-pla-kalman-fast-gpuopt-nodetrend.mq5`
  - Exibição exclusiva via `InpViewMode`: ou só ondas (FFT) ou só feed; buffers contrários são preenchidos com `EMPTY_VALUE`.
  - Removido o rótulo/feed label duplicado e qualquer traço do feed quando ondas estão visíveis.
  - Eliminado código de fallback CPU FFT; processamento permanece somente GPU (barra é ignorada se GPU indisponível).
  - `FeedTrace` redimensionado e preenchido apenas no modo feed; propriedades do plot ajustadas (único plot 9).
  - Feed ZigZag agora usa o timeframe selecionado (`InpFeedTimeframe`), convertendo a barra do gráfico para `shift` correto e invertendo a ordem para cronologia; evita paradas quando o gráfico está em outro timeframe (ex.: M15 com feed M1).
  - Pré-carregamento de histórico do feed: `HistorySelect` + `CopyClose` exigem pelo menos `InpFFTWindow` barras; se faltar histórico, loga WARN e não avança a barra até concluir o download.

## 2025-11-17
- `WaveSpecZZ_1.0.3-pla-kalman-fast-gpuopt-nodetrend.mq5`
  - Adicionado combo `InpFeedData` (PLA / ZigZag / Close) como fonte única de `feed_data`.
  - Removidos todos os fallbacks de feed: se PLA ou ZigZag falham, a barra é ignorada e loga erro explícito; modo inválido também aborta a barra.
  - Removidos fallbacks de tempo para ETA: se `PeriodSeconds` ou `delta` entre barras não forem válidos, ETA é desativada na barra e o erro é logado (sem usar 60s).
  - Mantidos logs claros indicando ausência de fallback.
