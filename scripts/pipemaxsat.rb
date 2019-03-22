#!/usr/bin/env ruby

#autcross with optional global timeout (the sanity checks can take long too)
tout = ARGV[0].to_i
ARGV.shift if tout > 0

tmpfile=`mktemp -p /dev/shm`.strip
["INT","TERM"].each{|sig| Signal.trap(sig) { `rm #{tmpfile}`; exit }}

cmd = "open-wbo -verbosity=0"
cmd = "timeout #{tout} "+cmd if tout>0

# puts cmd
# exit

def each_job(file, &fun)
  job = ""
  firstp = true

  file.each_line do |l|
    tok = l.split

    if tok[0] == "p"
      unless firstp
        fun.call(job)
        job = ""
      end
      firstp = false
    end

    job += l
  end
  #for the rest
  fun.call(job)
end

each_job(STDIN) do |a|
  open(tmpfile, 'w'){ |f| f.write a }
  # pipe = IO.popen("#{cmd} #{tmpfile}", 'r+')
  # pipe.close_write
  # STDOUT.write pipe.read
  puts `#{cmd} #{tmpfile}`
end

`rm #{tmpfile}`
