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
