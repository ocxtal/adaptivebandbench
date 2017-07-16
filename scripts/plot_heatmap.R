#! /usr/bin/env rscript

library(ggplot2)

plot_heatmap <- function(title, in_filename, out_filename) {
	data <- read.table(in_filename, head = T)
	label <- names(data)

	data[[label[1]]] <- as.character(data[[label[1]]])
	data[[label[2]]] <- as.character(data[[label[2]]])
	data$level <- cut(as.numeric(data[[label[3]]]),
		breaks = c(0, 90:101),
		labels = as.character(c(0, 90:100)),
		right = F)

	# color <- c(colorRampPalette(c("red", "blue"))(10), "white")

	# plot id_bw
	p <- ggplot(data, aes(
		x = data[[label[1]]],
		y = data[[label[2]]],
		value = data[[label[3]]],
		label = data[[label[3]]]))
	p <- p + geom_tile(aes(fill = data[[label[3]]]), colour = 'black') +
		geom_text(aes(label = data[[label[3]]])) +
		scale_fill_gradient(limits = c(80, 100), low = "gray", high = "white", na.value = "gray") +
		# scale_fill_gradient(limits = c(80, 100), low = "white", high = "blue", na.value = "red") +
		# scale_fill_manual(values = color, labels = as.character(c(0, 90:100))) +
		scale_x_discrete(limits = unique(data[[label[1]]]), expand = c(0, 0)) +
		scale_y_discrete(limits = unique(data[[label[2]]]), expand = c(0, 0)) +
		labs(x = label[1], y = label[2]) + 
		theme(axis.ticks = element_blank(), axis.line = element_line(size = 0.5), legend.position = "none") +
		ggtitle(title)
	ggsave(out_filename, p, width = 100, height = 80, units = 'mm')
}

plot_heatmap('Bandwidth-Identity', 'bw_id.txt', 'bw_id.eps')
plot_heatmap('Length-Identity (W = 16)', 'id_len_16.txt', 'id_len_16.eps')
plot_heatmap('Length-Identity (W = 32)', 'id_len_32.txt', 'id_len_32.eps')
plot_heatmap('X-Ge (Gi = 0)', 'x_ge_linear_2_x_y_2.txt', 'x_ge_linear_2_x_y_2.eps')
plot_heatmap('X-Gi (Ge = 1)', 'x_gi_affine_2_x_y_1.txt', 'x_gi_affine_2_x_y_1.eps')
plot_heatmap('X-Gi (Ge = 2)', 'x_gi_affine_2_x_y_2.txt', 'x_gi_affine_2_x_y_2.eps')
plot_heatmap('X-Gi (Ge = 3)', 'x_gi_affine_2_x_y_3.txt', 'x_gi_affine_2_x_y_3.eps')
plot_heatmap('Gi-Ge', 'gige_id_ddiag_extract.txt', 'gige_id.eps')

