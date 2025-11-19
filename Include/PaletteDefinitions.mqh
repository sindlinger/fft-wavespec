//+------------------------------------------------------------------+
//| PaletteDefinitions.mqh                                           |
//| Shared color palettes for WaveSpectrum indicators                |
//+------------------------------------------------------------------+
#ifndef __TFT_WAVESPECTRUM_PALETTE_DEFINITIONS_MQH__
#define __TFT_WAVESPECTRUM_PALETTE_DEFINITIONS_MQH__

enum COLOR_PRESET { PRESET_ELEGANT=0, PRESET_VIRIDIS=1, PRESET_PLASMA=2, PRESET_CIVIDIS=3, PRESET_SUNSET=4, PRESET_TOL=5, PRESET_MONO=6 };

inline color RGB(const int r, const int g, const int b)
{
   const uchar ur=(uchar)r;
   const uchar ug=(uchar)g;
   const uchar ub=(uchar)b;
   return (color)(((uint)ur) | (((uint)ug) << 8) | (((uint)ub) << 16));
}

const color PALETTE_VIRIDIS_VALUES[12] =
{
   RGB(68,1,84), RGB(71,44,122), RGB(59,81,139), RGB(44,113,142),
   RGB(33,144,141), RGB(39,173,129), RGB(92,200,99), RGB(150,219,64),
   RGB(208,226,36), RGB(244,229,38), RGB(254,231,51), RGB(241,229,103)
};

const color PALETTE_PLASMA_VALUES[12] =
{
   RGB(13,8,135), RGB(75,3,161), RGB(125,3,168), RGB(168,34,150),
   RGB(203,70,121), RGB(229,107,93), RGB(248,148,65), RGB(253,195,40),
   RGB(240,249,33), RGB(209,248,45), RGB(173,238,70), RGB(132,222,94)
};

const color PALETTE_CIVIDIS_VALUES[12] =
{
   RGB(0,32,76), RGB(0,48,113), RGB(0,63,133), RGB(53,81,134),
   RGB(95,99,132), RGB(136,119,127), RGB(175,142,120), RGB(208,168,108),
   RGB(233,198,93), RGB(247,229,81), RGB(249,242,144), RGB(236,245,191)
};

const color PALETTE_CARTO_SUNSET_VALUES[12] =
{
   RGB(4,58,74), RGB(32,89,103), RGB(67,120,127), RGB(107,147,146),
   RGB(152,174,159), RGB(192,190,162), RGB(224,184,153), RGB(244,165,143),
   RGB(244,129,122), RGB(232,91,104), RGB(202,52,103), RGB(160,26,99)
};

const color PALETTE_TOL_LIGHT_VALUES[12] =
{
   RGB(119,158,203), RGB(119,193,142), RGB(255,190,122), RGB(246,124,95),
   RGB(204,120,188), RGB(153,153,153), RGB(255,255,148), RGB(161,217,155),
   RGB(197,219,239), RGB(255,204,188), RGB(217,196,237), RGB(182,232,199)
};

struct SpectralMixDefinition
{
   double primary_nm;
   double secondary_nm;
   double primary_weight;
   double secondary_weight;
};

const SpectralMixDefinition g_spectral_palette_definitions[12] =
{
   {650.0,610.0,0.70,0.30}, {560.0,540.0,0.60,0.40}, {545.0,515.0,0.65,0.35},
   {498.0,470.0,0.60,0.40}, {575.0,555.0,0.60,0.40}, {650.0,440.0,0.55,0.45},
   {635.0,460.0,0.45,0.55}, {620.0,595.0,0.60,0.40}, {555.0,505.0,0.55,0.45},
   {508.0,486.0,0.50,0.50}, {590.0,570.0,0.55,0.45}, {470.0,450.0,0.65,0.35}
};

void GetPresetColors(COLOR_PRESET preset, color &out[])
{
   int i;
   switch(preset)
   {
      case PRESET_ELEGANT:
         out[0]=clrDarkSlateBlue; out[1]=clrSlateBlue; out[2]=clrRoyalBlue; out[3]=clrSteelBlue;
         out[4]=clrTeal; out[5]=clrDarkCyan; out[6]=clrSeaGreen; out[7]=clrMediumSeaGreen;
         out[8]=clrOliveDrab; out[9]=clrGoldenrod; out[10]=clrDarkOrange; out[11]=clrTomato;
         break;
      case PRESET_VIRIDIS: for(i=0;i<12;i++) out[i]=PALETTE_VIRIDIS_VALUES[i]; break;
      case PRESET_PLASMA:  for(i=0;i<12;i++) out[i]=PALETTE_PLASMA_VALUES[i];  break;
      case PRESET_CIVIDIS: for(i=0;i<12;i++) out[i]=PALETTE_CIVIDIS_VALUES[i]; break;
      case PRESET_SUNSET:  for(i=0;i<12;i++) out[i]=PALETTE_CARTO_SUNSET_VALUES[i]; break;
      case PRESET_TOL:     for(i=0;i<12;i++) out[i]=PALETTE_TOL_LIGHT_VALUES[i];  break;
      case PRESET_MONO:    for(i=0;i<12;i++){ int v=60+(i*10); out[i]=RGB(v,v,v); } break;
      default:             for(i=0;i<12;i++) out[i]=clrSilver; break;
   }
}

#endif // __TFT_WAVESPECTRUM_PALETTE_DEFINITIONS_MQH__
