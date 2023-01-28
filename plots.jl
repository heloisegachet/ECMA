using Plots

include("./data/10_ulysses_3.tsp")

partition = [[4,6,7,8], [5,9], [1,2,3,10]]

gr()
plot()
# draw points with clustering
for k in 1:K
	scatter!(coordinates[partition[k], 1], coordinates[partition[k], 2], series_annotations = text.(partition[k], :bottom), legend=false)
end
gui()