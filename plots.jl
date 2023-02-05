using Plots

include("./data/10_ulysses_3.tsp")

p_h = [[5, 9, 4], [10, 2, 1, 3], [7, 8, 6]]
val_h = 140.54478701953713
p_ls = [[5, 9], [10, 2, 1, 3], [7, 8, 6, 4]]
val_ls = 136.99527629589417

### solution before local search
no_ls = plot(show=true, fmt=:svg)
# draw points with clustering
for k in 1:K
	scatter!(coordinates[p_h[k], 1], coordinates[p_h[k], 2], series_annotations = text.(p_h[k], :bottom), legend=false)
end

savefig(no_ls, "sol_before_ls.svg")

## solution after local search
with_ls = plot(show=true, fmt=:svg)
# draw points with clustering
for k in 1:K
	scatter!(coordinates[p_ls[k], 1], coordinates[p_ls[k], 2], series_annotations = text.(p_ls[k], :bottom), legend=false)
end

savefig(with_ls, "sol_after_ls.svg")