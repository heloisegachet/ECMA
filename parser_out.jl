using CSV, XLSX
using DataFrames

function write(output_name, filename, time, sol, value, obj_lb, gap)
    filename_sol = string(output_name, "_sol.csv")
	filename_sol_xlsx = string(output_name, "_sol.xlsx")
	if isfile(filename_sol)
		df = CSV.read(filename_sol, DataFrame)
	else
		df = DataFrame("nom_fichier"=> [], 
						"time"=>[], 
						"value"=>[], "sol"=>[], "obj_bound"=>[],"gap"=>[])
		push!(df, [" "^30, 0, " "^30, " "^200, " "^30," "^30])
	end

    filename_info = split(filename, "/")[2]
	replace_row = false
	for row in eachrow(df)
		if row[:nom_fichier] == filename_info
			row[:time] = round(time, digits=3)
			value = split(string(value), " / ")
			if length(value) > 1
				first_val = round(parse(Float64, value[1]), sigdigits=8)
				second_val = round(parse(Float64, value[2]), sigdigits=8)
				value = string(first_val, " / ", second_val)
			elseif isa(value, Number)
				value = round(value[1], sigdigits=8)
			else
				value = value[1]
			end
			row[:value] = value
			row[:sol] = array_to_string(sol)
			row[:obj_bound] = string(obj_lb)
			row[:gap] = string(gap)
			replace_row = true
		end
	end
	if !replace_row
		value = split(string(value), " / ")
		if length(value) > 1
			first_val = round(parse(Float64,value[1]), sigdigits=8)
			second_val = round(parse(Float64,value[2]), sigdigits=8)
			value = string(first_val, " / ", second_val)
		elseif isa(value, Number)
			value = round(value[1], sigdigits=8)
		else
			value = value[1]
		end
		push!(df, [filename_info, time, string(value), array_to_string(sol), string(obj_lb),string(gap)])
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
	if length(result)>200
		return string(result[1:197], "...")
	else 
		return result
	end
end