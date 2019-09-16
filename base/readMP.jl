function readMP(filename; folder = "data", sep = ',')
   # 1 file
   if typeof(filename) != Array{String,2}
    df = CSV.read(joinpath(folder*"\\"*filename*".csv"); delim = sep)
    df = df[collect(x for x in names(df) if !ismissing(df[1,x]))] # Removing useless columns
    df = df[collect(x for x in 1:size(df,1) if !ismissing(df[x,end])),:] # Removing useless rows

    # Removing the ";" at the end
    rename!(df,names(df)[end] => Symbol(string(names(df)[end])[1:end-1]))
    for i in 1:size(df,1)
      df[i,end] = df[i,end][1:end-1]
    end
    return df
  else
    # More than 1 file
    DF = Dict()
    for i in 1:length(filename)
      df = CSV.read(joinpath(folder*"\\"*filename[i]*".csv"); delim = sep[i])
      df = df[collect(x for x in names(df) if !ismissing(df[1,x]))] # Removing useless columns
      df = df[collect(x for x in 1:size(df,1) if !ismissing(df[x,end])),:] # Removing useless rows

      # Removing the ";" at the end
      rename!(df,names(df)[end] => Symbol(string(names(df)[end])[1:end-1]))
      for j in 1:size(df,1)
        df[j,end] = df[j,end][1:end-1]
      end
      DF[i] = df
    end
    # return DF
    return [DF[i] for i in 1:length(filename)]
  end
end
