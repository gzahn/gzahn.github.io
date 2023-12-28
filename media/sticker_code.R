library(hexSticker)
library(ggplot2)
imgurl <- "./media/leaf_fungi.png"

s <- sticker(imgurl, package = "",#filename = "./media/bootcamp_hex.png",
        s_width = .8,
        s_height = .5,
        s_x = .9,
        s_y = 1,
        p_y = 1.5,
        p_size = 14,p_color = 'black',p_fontface = 'bold',p_family = 'arial',
        h_fill = "gray90",
        spotlight = TRUE,
        l_y = 1,
        l_alpha = .25,
        l_width = 100,
        l_height = 5)
s + annotate('text',x = 1,y=1.3,label="Boot Camp",size=16,fontface=2) +
    annotate('text',x = 1,y=1.5,label="Microbiome",size=16,fontface=2)
ggsave("./media/bootcamp_hex.png")
