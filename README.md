# 🎨 ProteoFigure Studio - Publication-Grade Figure Maker

**ProteoFigure Studio**, diferansiyel proteomik analizleri için tasarlanmış yüksek kaliteli bir görselleştirme ve yayın figürü üreticisidir.

Bu modül, serideki üç ana pipeline'ın (**Proteomics_Pipeline.Rmd**, **Ultimate_Functional_Analysis.Rmd** ve **Proteomic_Footprinting.Rmd**) ürettiği tüm veri çıktılarını doğrudan içe aktararak yalnızca klasik grafiklerle sınırlı kalmadan; zengin grafik seçenekleri, prestijli dergi renk paletleri ve tam manuel özelleştirme imkanı sunar.

---

## 🌟 Neden ProteoFigure Studio?

- **Sadece Dot Plot Değil:** Fonksiyonel zenginleştirme için **Dot Plot**, **Yatay/Dikey Barplot**, **Lollipop Chart**, **GSEA Şelale Grafiği**.
- **Diferansiyel Ekspresyon için Zengin Seçenekler:** **Enhanced Volcano Plot** (`ggrepel` ile arama yaparak istenen proteinleri etiketleme), **Tekil Protein Boxplot / Violin**, **PCA Güven Elips Dağılımı**.
- **Transkripsiyon Faktörleri (TF):** Çift yönlü **Diverging Barplot**, **TF Lollipop** ve skor dağılımları.
- **Prestijli Renk Paletleri:** Nature Publishing Group (NPG), NEJM, The Lancet, JAMA, JCO, Science (AAAS), Viridis, ColorBrewer ve Özel Renk Seçiciler (Hex).
- **Yayın Formatlarında Kaydetme:** Vektörel **PDF**, **SVG** (Illustrator uyumlu), **300 & 600 DPI PNG**, **TIFF** formatlarında tek tıkla indirme; Nature tek sütun (89 mm) ve çift sütun (183 mm) boyut önayarları.

---

## 🚀 Nasıl Kullanılır?

### 1. Yöntem: İnteraktif Arayüz (Shiny App) - Tavsiye Edilen 🖥️

Uygulamayı başlatmak için RStudio veya R konsolunda aşağıdaki komutu çalıştırmanız yeterlidir:

```r
source("c:/Users/serife-ayaz/Desktop/Differential_Proteomics_Pipeline/Figure_Maker/run_figure_maker.R")
```

Veya `In_Silico_Protocols` projesindeyseniz:
```r
source("run_figure_maker.R")
```

Açılan modern arayüzde:
1. Sonuç klasörünüzü (`OMV_Treatment_Sec` vb.) seçin veya dosyanızı yükleyin.
2. Grafik türünü, metrikleri (GeneRatio, Count, FoldEnrichment) ve paleti belirleyin.
3. Canlı önizlemede grafiğinizi inceleyin.
4. "Download Figure" butonuna basarak dergi formatında çıktınızı alın!

### 2. Yöntem: RMarkdown / Script İçinde Kodla Kullanım 📝

Arayüz kullanmak istemiyorsanız, doğrudan `figure_maker_engine.R` fonksiyonlarını projenize çağırabilirsiniz:

```r
source("figure_maker_engine.R")

# 1. Zenginleştirme Dot Plot
df_go <- read.csv("Ultimate_GO_BP_Enriched_Young_OMV_vs_Young_Control.csv")
p_dot <- plot_enrichment_dot(df_go, x_metric = "GeneRatio", color_metric = "p.adjust", top_n = 15)

# 2. Volcano Plot (Özel Gen Etiketleriyle)
df_dea <- read.csv("all_proteins_Young_OMV_vs_Young_Control.csv")
p_volc <- plot_volcano(df_dea, selected_genes = c("IL6", "C3", "CASP3"), label_type = "custom")

# 3. Yüksek Kalitede PDF Kaydet
export_figure(p_dot, file_path = "Figure1.pdf", format = "pdf", width = 7.2, height = 6.0)
```

Detaylı örnek kodlar ve çıktılar için [Figure_Maker_Tutorial.Rmd](Figure_Maker_Tutorial.Rmd) dosyasına göz atabilirsiniz.

---

## 📁 Dosya Yapısı

```text
Figure_Maker/
├── app.R                     # İnteraktif Shiny App
├── figure_maker_engine.R     # Çekirdek modüler grafik fonksiyonları motoru
├── run_figure_maker.R        # Tek tıklamayla başlatıcı script
├── Figure_Maker_Tutorial.Rmd # Kodla çalışan RMarkdown kılavuzu
└── README.md                 # Bu dokümantasyon
```
