#!/usr/bin/env ruby
#filter for automata that determinize quickly and save slow in other file
require 'open3'
tout = ARGV[0].to_i
slowdet_out = ARGV[1]
cmd = "timeout #{tout} autfilt -D -P --high"

def each_aut(file, &fun)
  aut = ""
  file.each_line do |l|
    aut += l
    l.strip!
    aut = "" if l == "--ABORT--"
    if l == "--END--"
      fun.call(aut)
      aut = ""
    end
  end
end

each_aut(STDIN) do |a|
  Open3.popen3(cmd) do |stdin, stdout, stderr, wait_thr|
    stdin.puts a
    stdin.close
    stdout.read
    stderr.read

    exit_status = wait_thr.value.exitstatus

    if exit_status != 124 #no timeout
      puts a
    else
      open(slowdet_out, 'a'){ |f| f.puts(a) }
    end

  end
end

