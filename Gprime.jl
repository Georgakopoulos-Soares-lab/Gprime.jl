# Copyright Congzhou M Sha (cms6712@psu.edu)
# All Rights Reserved

using Pkg
Pkg.activate(".")
using ArgParse
using DataFrames
using CSV
using Serialization
using ProgressBars
using GZip
using Base.Threads
using JLD2

println("Using $(nthreads()) threads")

const alphabet = Set(['A', 'G', 'C', 'T'])
const k=12

@inline @fastmath function translate_RNA(nt)
    if nt == 'A'
        return '0'
    elseif nt == 'G'
        return '1'
    elseif nt == 'C'
        return '2'
    else
        return '3'
    end
end

@inline @fastmath function base_to_int(nt)
    if nt == 'A'
        return 0
    elseif nt == 'G'
        return 1
    elseif nt == 'C'
        return 2
    else
        return 3
    end
end

@inline @fastmath function complement(nt)
    if nt == '0'
        return '3'
    elseif nt == '1'
        return '2'
    elseif nt == '2'
        return '1'
    else
        return '0'
    end
end

@inline @fastmath function bit_encode(str, k=k)
    base_4_rep = join(translate_RNA.(collect(str)));
    parse(UInt64, base_4_rep; base=4)
end

@inline @fastmath function bit_encode_complement(str, k=k)
    base_4_rep = join(complement.(translate_RNA.(collect(reverse(str)))));
    parse(UInt64, base_4_rep; base=4)
end

@inline @fastmath function get_kmers(s, k=k)
    Set(bit_encode(s[i:i+k-1]) for i in 1:length(s)-k+1 if 
                all(j in alphabet for j in s[i:i+k-1]))
end

@inline @fastmath function get_set(in_f, k=k, SET_TYPE=BitSet)
    mask = 1 << (k << 1) - 1
    result = SET_TYPE()
    sizehint!(result, 1 << (k << 1))
    prev = ""
    while !eof(in_f)
        line = readline(in_f);
        @inbounds prev_char = line[1];
        for i in 1:200
            if eof(in_f) || prev_char != '>'
                break
            end
            next_line = readline(in_f)
            prev_char = next_line[1]
            line = line * next_line
        end
        if length(line) == 0 || line[1] == '>'
            prev = "";
            continue
        end
        prev = prev * line;
        if length(prev) < k
            continue
        end
        
        translated = bit_encode(prev[1:k])
        for i in (k+1):length(prev)
            @inbounds translated = ((translated << 2) | base_to_int(prev[i])) & mask
            push!(result, translated)
        end
        @inbounds prev = prev[length(prev)-k+2:end]
    end
    return result
end

@inline @fastmath function get_sets(genome_path="genomes", out_path=".", k=k, SET_TYPE=BitSet, mem_reduce=false)
    fns = [i for i in readdir(genome_path) if occursin("fasta", i)]
    println(fns)
    fns = [i for i in fns if i[end-4:end] == "fasta" || i[end-7:end] == "fasta.gz"]
    Base.Filesystem.mkpath("$(out_path)/genome_$(k)mers");
    Threads.@sync for fn in fns
        Threads.@spawn begin
            in_fn = "$(genome_path)/$(fn)"
            stub = split(fn, "fasta")[1]
            out_fn = "$(out_path)/genome_$(k)mers/$(stub)ser";
            if !Base.Filesystem.ispath(out_fn)
                println("Processing $(in_fn)")
                in_f = in_fn[end-1:end] == "gz" ? GZip.open(in_fn) : open(in_fn)
                if !mem_reduce
                    in_f = IOBuffer(read(in_f))
                end
                result = get_set(in_f, k, SET_TYPE)
                close(in_f)

                if SET_TYPE == BitSet && length(result) < ((1 << (k << 1))) / 25
                    temp = Set(result)
                    if Base.summarysize(temp) < Base.summarysize(result)
                        result = Set(result)
                    end
                end
            
                open(out_fn, "w") do f
                    serialize(out_fn, result)
                end
            end
        end
    end
end

@inline @fastmath function get_primes(out_path=".", k=k)
    path = "$(out_path)/genome_$(k)mers"
    prime_path = "$(out_path)/quasi_$(k)_primes"
    Base.Filesystem.mkpath(prime_path)
    fns = [i for i in readdir(path) if occursin("ser", i)]
    fns = [i for i in fns if i[end-2:end] == "ser"]
    genomes = Array{Any}(undef, length(fns))
    @sync for (i, fn) in enumerate(fns)
        @spawn begin
            println("Loading $(fn)...")
            open("$(path)/$(fn)") do f
                genomes[i] = (fn[1:end-3], deserialize(f))
            end
        end
    end
    
    if nthreads() > 1
        Threads.@sync for (k1, v1) in genomes
            Threads.@spawn begin
                write_path = "$(prime_path)/p_$(k1)ser"
                if !Base.Filesystem.ispath(write_path)
                    result = copy(v1)
                    println("Isolating $(k1)")
                    println(write_path)
                    for (k2, v2) in genomes
                        if k1 == k2
                            continue
                        end
                        setdiff!(result, v2)
                    end
                    println("Saving $(write_path)...")
                    open(write_path, "w") do f
                        serialize(f, result)
                    end
                end
            end
        end
    else
        for (k1, v1) in ProgressBar(genomes)
            write_path = "$(prime_path)/p_$(k1)ser"
            if !Base.Filesystem.ispath(write_path)
                result = copy(v1)
                println("Isolating $(k1)")
                for (k2, v2) in genomes
                    if k1 == k2
                        continue
                    end
                    result = setdiff(result, v2)
                end
                open(write_path, "w") do f
                    serialize(f, result)
                end
            end
        end
    end
end

@fastmath function get_primes_mem_reduce(out_path=".", k=k)
    path = "$(out_path)/genome_$(k)mers"
    prime_path ="$(out_path)/quasi_$(k)_primes"
    Base.Filesystem.mkpath(prime_path)

    fns = [i for i in readdir(path) if occursin("ser", i)]
    fns = [i for i in fns if i[end-2:end] == "ser"]

    if nthreads() > 1
        Threads.@sync for i in 1:length(fns)
            Threads.@spawn begin
                fn = fns[i]
                k1 = fn[1:end-3]
                write_path = "$(prime_path)/p_$(k1)ser"
                if !Base.Filesystem.ispath(write_path)
                    println("Isolating $(k1)")
                    result = nothing
                    open("$(path)/$(fn)") do f
                        result = deserialize(f)
                    end
           
                    for fn2 in fns
                        if fn2 == fn
                            continue
                        end
                        k2 = fn2[1:end-3]
                        println("Loading $(path)/$(fn2)")
                        v2 = nothing
                        open("$(path)/$(fn2)") do f
                            v2 = deserialize(f)
                        end
                        setdiff!(result, v2)
                    end
                    println("Saving $(write_path)...")
                    open(write_path, "w") do f
                        serialize(f, result)
                    end
                end
            end
        end
    else
        for fn in ProgressBar(fns)
            k1 = fn[1:end-3]
            write_path = "$(prime_path)/p_$(k1)ser"
            if !Base.Filesystem.ispath(write_path)
                println("Isolating $(k1)")
                result = nothing
                open("$(path)/$(fn)") do f
                    result = deserialize(f)
                end
                for fn2 in fns
                    if fn2 == fn
                        continue
                    end
                    k2 = fn2[1:end-3]
                    v2 = nothing
                    open("$(path)/$(fn2)") do f
                        v2 = deserialize(f)
                    end
                    setdiff!(result, v2)
                end
                open(write_path, "w") do f
                    serialize(f, result)
                end
            end
        end
    end
end

@inline @fastmath function create_master_table(out_path=".", k=k)
    master_path = "$(out_path)/master_table_$(k).ser"
    if Base.Filesystem.ispath(master_path)
        tab = nothing
        open(master_path) do f
            tab = deserialize(f)
        end
        species = tab["species"]
        master_table = tab["master_table"]
        return species, master_table
    end
    println("Creating master table...")
    prime_path = "$(out_path)/quasi_$(k)_primes"
    primes = [fn for fn in readdir(prime_path) if fn[end-2:end] == "ser"]
    species = sort([i[3:end-4] for i in primes])
    N = length(species) + 1
    INDEX_TYPE = N <= 256 ? UInt8 : (N <= 65536 ? UInt16 : (N <= (1 << 32) ? UInt32 : UInt64))
    rev_species = Dict(j => INDEX_TYPE(i) for (i,j) in enumerate(species))
    master_table = zeros(INDEX_TYPE, 1 << (k << 1))
    for fn in primes
        s = fn[3:end-4]
        species_id = rev_species[s]
        open("$(prime_path)/$(fn)") do f
            for p in deserialize(f)
                master_table[p+1] = species_id
            end
        end
    end
    open(master_path, "w") do f
        serialize(f, Dict("species" => species, "master_table" => master_table))
    end
    return species, master_table
end

@fastmath function classify(in_f, master_table, k=k)
    mask = 1 << (k << 1) - 1
    if nthreads() > 1
        read_counts = zeros(UInt64, 256)
        prev = ""
        while !eof(in_f)
            line = readline(in_f);
            @inbounds prev_char = line[1];
            for i in 1:200
                if eof(in_f) || prev_char == '>'
                    break
                end
                next_line = readline(in_f)
                if length(next_line) == 0
                    break
                end
                prev_char = next_line[1]
                line = line * next_line
            end
            if length(line) == 0 || line[1] == '>'
                prev = "";
                continue
            end
            prev = prev * line;
            if length(prev) < k
                continue
            end
            
            translated = bit_encode(prev[1:k])
            for i in (k+1):length(prev)
                @inbounds translated = ((translated << 2) | base_to_int(prev[i])) & mask
                read_counts[master_table[translated+1]+1] += 1
            end
            @inbounds prev = prev[length(prev)-k+2:end]
        end
        return read_counts
    end
    
    pbar = ProgressBar()
    read_counts = zeros(UInt64, 256)
    prev = ""
    while !eof(in_f)
        line = readline(in_f)
        @inbounds prev_char = line[1];
        for i in 1:200
            if eof(in_f) || prev_char != '>'
                break
            end
            next_line = readline(in_f)
            prev_char = next_line[1]
            line = line * next_line
        end
        if length(line) == 0 || line[1] == '>'
            prev = "";
            continue
        end
        prev = prev * line;
        if length(prev) < k
            continue
        end

        update(pbar)
        translated = bit_encode(prev[1:k])
        for i in (k+1):length(prev)
            @inbounds translated = ((translated << 2) | base_to_int(prev[i])) & mask
            read_counts[master_table[translated+1]+1] += 1
        end
        @inbounds prev = prev[length(prev)-k+2:end]
    end
    return read_counts
end

@inline @fastmath function classify_all(master_table, sample_path="samples", out_path=".", k=k)
    Base.Filesystem.mkpath("$(out_path)/classified");
    fns = [i for i in readdir(sample_path) if occursin("fasta", i)]
    fns = [i for i in fns if i[end-4:end] == "fasta" || i[end-7:end] == "fasta.gz"]
    Threads.@sync for fn in fns
        Threads.@spawn begin
            in_fn = "$(sample_path)/$(fn)"
            println("Classifying $(in_fn)...")
            stub = split(fn, "fasta")[1]
            in_f = in_fn[end-1:end] == "gz" ? GZip.open(in_fn) : open(in_fn)
            result = classify(in_f, master_table, k)
            close(in_f)
            println("Saving $(out_path)/classified/$(stub)h5")
            jldopen("$(out_path)/classified/$(stub)h5", "w") do f
                f["counts"] = result
            end
        end
    end
end

function parse_commandline()
    s = ArgParseSettings("This program implements the entire pipeline from creation of a master table of the quasi-primes for a set of genomes to classifying sequencing data.")

    @add_arg_table s begin
        "-k"
            help = "the length of kmers you wish to classify"
            arg_type = Int
            default = 12
        "--genome_path", "-g"
            help = "the folder containing the reference genomes with respect to which the master table will be created (*.fasta or *.fasta.gz)"
            arg_type = String
            default = "genomes"
        "--sample_path", "-s"
            help = "the folder containing the sequencing data you wish to classify (*.fasta or *.fasta.gz)"
            arg_type = String
            default = "samples"
        "--out_path", "-o"
            help = "the folder where all intermediate data will be stored"
            arg_type = String
            default = "."
        "--mem_reduce"
            help = "reduce memory usage at the cost of increased disk reading while computing the primes. Can be used independently of the --bitsets flag."
            action = :store_true
        "--clean", "-c"
            help = "a flag which indicates if you wish to clean the current directory of intermediate files. Otherwise, the program resumes its progress."
            action = :store_true
        "--bitsets", "-b"
            help = "a flag which indicates if you wish to use BitSets vs Sets in Julia. Note that BitSets are more CPU-efficient for dense data, such as if you choose a length k such that most sequences will be assigned to a taxon; Sets are more memory- and CPU-efficient for sparse data."
            action = :store_true
        "--skip_sets"
            help = "skip computing the sets of kmers"
            action = :store_true
        "--skip_primes"
            help = "skip computing the primes"
            action = :store_true
        "--skip_master"
            help = "skip computing the master table"
            action = :store_true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    clean = parsed_args["clean"]
    k = parsed_args["k"]
    genome_path = parsed_args["genome_path"]
    sample_path = parsed_args["sample_path"]
    out_path = parsed_args["out_path"]
    use_bitsets = parsed_args["bitsets"]
    mem_reduce = parsed_args["mem_reduce"]
    skip_sets = parsed_args["skip_sets"]
    skip_primes = parsed_args["skip_primes"]
    skip_master = parsed_args["skip_master"]
    SET_TYPE = use_bitsets ? BitSet : Set{UInt64}
    if clean
        rm("$(out_path)/genome_$(k)mers/*ser", force=true, recursive=true)
        rm("$(out_path)/quasi_$(k)_primes/*ser", force=true)
        rm("$(out_path)/master_table_$(k).ser", force=true)
    end
    if !skip_sets
        get_sets(genome_path, out_path, k, SET_TYPE, mem_reduce)
    end
    
    if !skip_primes
        if mem_reduce
            get_primes_mem_reduce(out_path, k)
        else
            get_primes(out_path, k)
        end
    end
    species, master_table = create_master_table(out_path, k)
    jldopen("species.h5", "w") do f
        f["species"] = species
    end
    classify_all(master_table, sample_path, out_path, k)
end

main()
