library(hexSticker)
library(tidyverse)
library(sysfonts)

# create dag to use in sticker ------------------------------------------------

dag <- ggplot() +
  # create nodes
  annotate("text", x = 1, y = 1, label= "X",
           col = "white", size = 3, fontface="bold") +
  annotate("text", x = 3, y = 1, label = "Y",
           col = "white", size = 3, fontface="bold") +
  annotate("text", x = 2, y = 2, label = "MEDIATOR",
           col = "skyblue", size = 5, fontface="bold") +
  # create edges
  ## x -> y
  geom_segment(aes(x = 1.2, y = 1, xend = 2.8, yend = 1),
               arrow = arrow(length = unit(0.03, "npc")),
               col = "white", size = 0.5, lineend="round") +
  ## x -> mediator
  geom_segment(aes(x = 1, y = 1.2, xend = 1.9, yend = 1.8),
               arrow = arrow(length = unit(0.03, "npc")),
               col = "white", size = 0.5, lineend="round") +
  ## mediator -> y
  geom_segment(aes(x = 2.1, y = 1.8, xend = 2.9, yend = 1.2),
               arrow = arrow(length = unit(0.03, "npc")),
               col = "white", size = 0.5, lineend="round") +
  # limit plot size
  xlim(c(0.5, 3.5)) + ylim(c(0.5, 2.5)) +
  # adjust background
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "grey25")) + theme_transparent()

# create hex sticker ----------------------------------------------------------

hexSticker::sticker(dag, s_width = 1.7, s_height = 1.7,  s_x = 1, s_y = 1,
                    h_fill = "grey25", package = "", h_color = "skyblue",
                    url = "www.gerkelab.com", u_color = "white", u_size = 1.3,
                    filename="man/figures/hex.png")
