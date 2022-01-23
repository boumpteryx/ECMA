# parser

function read_instance(MyFileName::String)
  path = "./Instances_ECMA/" * MyFileName
	# Si le fichier path existe
	if isfile(path) # "./Instances_ECMA/20_USA-road-d.BAY.gr"
		# L’ouvrir
		myFile = open(path) # "./Instances_ECMA/20_USA-road-d.BAY.gr"

    # Lire les premieres lignes
    n = parse(Int64, split(readline(myFile), " ")[3])
    s = parse(Int64, split(readline(myFile), " ")[3])
    t = parse(Int64, split(readline(myFile), " ")[3])
    S = parse(Int64, split(readline(myFile), " ")[3])
    d1 = parse(Int64, split(readline(myFile), " ")[3])
    d2 = parse(Int64, split(readline(myFile), " ")[3])

    p = Vector{Int64}(undef,0)
    line = split(readline(myFile), " ")
    number = line[3]
    number = chop(number, head=1,tail=1)
    number = parse(Int64, number)
    append!(p, number)
    popfirst!(line)
    popfirst!(line)
    popfirst!(line)
    for chunk in line
      chunk = chop(chunk)
      append!(p, parse(Int64, chunk))
    end

    ph = Vector{Int64}(undef,0)
    line = split(readline(myFile), " ")
    number = line[3]
    number = chop(number, head=1,tail=1)
    number = parse(Int64, number)
    append!(ph, number)
    popfirst!(line)
    popfirst!(line)
    popfirst!(line)
    for chunk in line
      chunk = chop(chunk)
      append!(ph, parse(Int64, chunk))
    end

    readline(myFile)
    # valeur 0 si il n'y a pas d'arc
    Mat_d = Array{Int64,2}(zeros(n,n))
    Mat_D = Array{Float64,2}(zeros(n,n))
		# Lire toutes les lignes d’un fichier
		data = readlines(myFile)
		for datum in data
      line = split(datum, " ")
      # pop!(line[3])
      line[3] = chop(line[3])
      # pop!(line[4])
      line[4] = chop(line[4])
      Mat_d[parse(Int64,line[1]),parse(Int64,line[2])] = parse(Int64,line[3])
      Mat_D[parse(Int64,line[1]),parse(Int64,line[2])] = parse(Float64,line[4])
		end

		return n, s, t, S, d1, d2, p, ph, Mat_d, Mat_D
	end
end
