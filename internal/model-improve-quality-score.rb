## for longer homopolymer runs, make the first and last bases be the lowest quality


fa = ARGV[0]
qual = ARGV[1]

$order= []
$reads = {}
name = ''

File.new(fa, 'r').each do |line|
  if line=~ /^>(\S+)/
    name = $1
    $reads[name] = {}
    $order << name
    
    $reads[name][:desc] = line.chomp
    $reads[name][:seq] = ''
    $reads[name][:qual] = []
  else
    $reads[name][:seq] << line.chomp
  end
end


File.new(qual, 'r').each do |line|
  if line=~ /^>(\S+)/
    name = $1
  else
    $reads[name][:qual] << line.chomp.split(/\s+/)
  end
end

# assign the qual value of the first base to the middle base, then assign the average qual of the last two bases to both ends of the homopolymer

$order.each do |r|
  desc = $reads[r][:desc]
  seq = $reads[r][:seq].split('')
  qualarray = $reads[r][:qual]

  qualqueue = []
  homoqueue = []
  homoqueue << seq.shift
  qualqueue = qualarray.shift

  puts ">#{desc}"
  0.upto(seq.size - 1) do |i|
    if seq[i] == homoqueue[-1]
      homoqueue << seq[i]
      qualqueue << qualarray[i]
    else # the end of a homoqueue
      headqual = qualqueue[0]
      if homoqueue.size == 1
        print "#{qualqueue[0]} "
      elsif homoqueue.size == 2
        print "#{qualqueue[0]} #{qualqueue[1]} "
      elsif homoqueue.size == 3 
        q1 = qualqueue[1]
        q2 = qualqueue[0]
        print "#{q1} #{q2} #{qualqueue[2]} "
      elsif homoqueue.size > 3
        qmid = qualqueue[0]
        qlasttwoavg = (qualqueue[-1] + qualqueue[-2])/2
        midp = (qualqueue.size -1) / 2  # the middle point

        shift = (qualqueue[midp] - qualqueue[-2])/ (qualqueue.size - 3)

        qualqueue[midp] = qmid
        qualqueue[0] = qlasttwoavg
        qualqueue[-1] = qlasttwoavg
        1.upto(qualqueue.size - 3) do |i|
          if i!= midp
            qualqueue[i] += shift
          end
        end

        print qualqueue.join(" ") + " "

        qualqueue = []
        homoqueue = []
        homoqueue << seq[i]
        qualqueue << qualarray[i]

      end
    end
  end
  print "\n"
end
        
        

      


