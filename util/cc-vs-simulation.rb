## 

realcc = ARGV[0]

simufof = ARGV[1]

hash = {}
avgs = {}

File.new(simufof, 'r').each do |line|
  if line=~/^(\S+)/
    ff = $1
    File.new(ff,'r').each do |j|
      if j=~ /^(\S+)\s+(\S+)\s+(\d+)/
        gene, cc, ss = $1, $2, $3.to_i

        if cc != 'NaN'
          hash[gene] = []  unless hash.key?(gene)
          hash[gene] << cc.to_f
        end
      end
    end
  end
end

hash.each_key do |gene|
  t = 0
  maxcc = 0
  hash[gene].each do |j|
    t += j
    if j.abs > maxcc
      maxcc = j
    end
  end
  avgs[gene] = [t / hash[gene].size.to_f, maxcc]
end


File.new(realcc, 'r').each do |line|
  if line=~ /^(\S+)\s+(\S+)\s+(\d+)/
    gene, cc, ss = $1, $2, $3.to_i
    if avgs.key?(gene) and cc!= 'NaN'
      puts "#{gene}\t#{avgs[gene][0]}\t#{cc}\t#{avgs[gene][1]}"
    end
  end
end
