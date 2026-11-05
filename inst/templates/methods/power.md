<!-- refs: edger, hart -->
<!-- packages: edgeR -->
The biological coefficient of variation was estimated with edgeR on the filtered genes (common BCV with a bootstrap over genes, {{p.bootstrap_replicates}} replicates). Indicative sample sizes per group for a two-group test on negative-binomial counts were computed with the approximation of Hart et al. over a grid of coefficients of variation ({{p.grid.cv}}), mean counts ({{p.grid.mu}}), fold changes ({{p.grid.fold_change}}) and significance levels ({{p.grid.alpha}}) at power {{p.grid.power}}.
