using CSV, XLSX
using DataFrames

function write(output_name, filename, time, sol, value, obj_lb, gap)
    filename_sol = string(output_name, "_sol.csv")
	filename_sol_xlsx = string(output_name, "_sol.xlsx")
	if isfile(filename_sol)
		df = CSV.read(filename_sol, DataFrame)
	else
		df = DataFrame("nom_fichier"=> [], "time"=>[], "value"=>[], "sol"=>[], "obj_bound"=>[],"gap"=>[])
	end

    filename_info = split(filename, "/")[2]
	replace_row = false
	for row in eachrow(df)
		if row[:nom_fichier] == filename_info
			row[:time] = time
			row[:value] = value
			row[:sol] = array_to_string(sol)
			row[:obj_bound] = obj_lb
			row[:gap] = gap
			replace_row = true
		end
	end
	if !replace_row
		push!(df, [filename_info, time, value, array_to_string(sol), obj_lb,gap])
	end

	CSV.write(filename_sol, df)
	XLSX.openxlsx(filename_sol_xlsx, mode="w") do xf
		sheet = xf[1]
		XLSX.writetable!(sheet, collect(DataFrames.eachcol(df)), DataFrames.names(df))
	end
end


function array_to_string(array)
	result = ""
	for i in 1:size(array)[1]
		result = string(result, "[")
		for j in 1:length(array[i])
			if (j < length(array[i]))
				result = string(result, array[i][j], ", ")
			else
				result = string(result, array[i][j])
			end
		end
		result = string(result, "]")
	end
	if length(result)>30
		return string(result[1:27], "...")
	else 
		return result
	end
end