function save_energy(energy::Array)
    dir  = "data"
    filename = "energy.dat"

    file = open(joinpath(dir,filename),"w")
    write(file, "# ENERGY for a ∈ [$a_begin,$a_end], c ∈ [$c_begin,$c_end] \n")

    for i in eachindex(a)
        for j in eachindex(c)
            write(file,"$(energy[i,j])\t")
        end
        write(file,"\n")
    end

    close(file)
end

function save_variance(variance::Array)
    dir  = "data"
    filename = "variance.dat"

    file = open(joinpath(dir,filename),"w")
    write(file, "# variance for a ∈ [$a_begin,$a_end], c ∈ [$c_begin,$c_end] \n")

    for i in eachindex(a)
        for j in eachindex(c)
            write(file,"$(variance[i,j])\t")
        end
        write(file,"\n")
    end

    close(file)
end

function load_data(filename::String)
    lines = readlines(filename)
    y_size = length(lines)
    x_size = 0
    i = 1
    
    while(true)
        if(lines[1][begin] == '#')
            continue
        end
        x_size = length(lines[i])
        i+=1
        break
    end
end
