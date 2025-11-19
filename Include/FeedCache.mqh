#ifndef __WAVESPEC_FEEDCACHE_MQH__
#define __WAVESPEC_FEEDCACHE_MQH__

// FeedCache.mqh — cache incremental de barras (Close) para qualquer indicador/EA
// Uso básico:
//   1) #include "Include\\FeedCache.mqh"
//   2) Declare FeedCache cache; int delta=0; bool from_file=false;
//   3) if(!EnsureFeedCache(cache, _Symbol, PERIOD_M1, need, true, "MeuPrefixo", delta, from_file)) return;
//      - Se já existir cache suficiente, não baixa nada.
//      - Se faltar barras, baixa só o delta (em chunks) e anexa, salvando de volta.
//      - Se símbolo/TF mudou ou não há cache, baixa tudo.
//   4) Use cache.close[idx] com ArraySetAsSeries=true (mais recente em [0]).
// Parâmetros:
//   enable_cache: true grava/ lê arquivo em FILE_COMMON; false só mantém em memória.
//   prefix: string para nomear o arquivo (evita conflito entre vários indicadores).
// Saídas:
//   delta_added: barras novas baixadas.
//   from_file_out: true se usou dados carregados de arquivo (mesmo com append).

struct FeedCache
{
    string          symbol;
    ENUM_TIMEFRAMES tf;
    double          close[];
    bool            loaded;
    bool            from_file;
};

// prefix: identifica o produto (ex.: "WaveSpecZZ", "MyEA"); evita conflito entre indicadores
string FeedCacheFileName(const string prefix, const string symbol, const ENUM_TIMEFRAMES tf)
{
    return StringFormat("%s_cache_%s_%s.bin", prefix, symbol, EnumToString(tf));
}

// Garante cache incremental (append em chunks) e salva se habilitado.
bool EnsureFeedCache(FeedCache &cache,
                     const string symbol,
                     const ENUM_TIMEFRAMES tf,
                     const int needed_bars,
                     const bool enable_cache,
                     const string prefix,
                     int &delta_added,
                     bool &from_file_out)
{
    delta_added=0;
    from_file_out=false;

    // 1) tenta carregar de arquivo se ainda não carregado
    if(enable_cache && !cache.loaded)
    {
        int h = FileOpen(FeedCacheFileName(prefix, symbol, tf), FILE_READ|FILE_BIN|FILE_SHARE_READ|FILE_COMMON);
        if(h!=INVALID_HANDLE)
        {
            int cnt = (int)FileReadInteger(h, INT_VALUE);
            if(cnt>0)
            {
                ArrayResize(cache.close, cnt);
                FileReadArray(h, cache.close, 0, cnt);
                cache.symbol = symbol;
                cache.tf     = tf;
                cache.loaded = true;
                cache.from_file = true;
                from_file_out = true;
            }
            FileClose(h);
        }
    }

    ArraySetAsSeries(cache.close, true);
    int cached = ArraySize(cache.close);
    bool same_feed = (cache.symbol==symbol && cache.tf==tf);

    if(!same_feed)
    {
        ArrayResize(cache.close, 0);
        cached = 0;
    }

    // 2) append em chunks até atingir needed_bars
    int max_chunk = 100000;
    double tmp[];
    ArraySetAsSeries(tmp, true);
    while(cached < needed_bars)
    {
        int want = MathMin(max_chunk, needed_bars - cached);
        ArrayResize(tmp, want);
        int got = CopyClose(symbol, tf, cached, want, tmp);
        if(got <= 0)
            break;
        int new_size = cached + got;
        ArrayResize(cache.close, new_size);
        for(int i=0;i<got;i++)
            cache.close[cached+i] = tmp[i];
        cached += got;
        delta_added += got;
    }

    cache.symbol = symbol;
    cache.tf     = tf;
    cache.loaded = (cached > 0);

    if(enable_cache && cache.loaded)
    {
        int h = FileOpen(FeedCacheFileName(prefix, symbol, tf), FILE_WRITE|FILE_BIN|FILE_COMMON);
        if(h!=INVALID_HANDLE)
        {
            FileWriteInteger(h, ArraySize(cache.close), INT_VALUE);
            FileWriteArray(h, cache.close, 0, ArraySize(cache.close));
            FileClose(h);
        }
    }

    // true somente se alcançamos o necessário
    return (cached >= needed_bars);
}

#endif // __WAVESPEC_FEEDCACHE_MQH__
