#!/usr/bin/env ruby
# -*- coding: utf-8 -*-

#
# lut_plot.rb
#
# Author: Rintaro Okamura
#
# Description:
#   This script makes a plot of look-up table.
#
# Dependency:
#   NArray gem
#   gnuplot gem
#
# Changelog:
#   20151118: First version
#   20151128: Add color and labels
#   20151129: Revise line styles and save filetype
#   20160117: Add surface albedo
#

require 'narray'
require 'gnuplot'
require 'open3'
require 'optparse'

MICORE_BIN = '../../src/micore'

params = ARGV.getopts('', 'iterplot:')

if ARGV.size < 5 then
  STDERR.puts <<heredoc
Usage: ruby #{$0} [--iterplot] input output albedo ref1 ref2
  input:  lut filename
  output: output png filename
  albedo: surface albedo
  ref1, ref2: reflectances
  [--iterplot]: make a plot of tau, cder, cost, ref1, ref2 for each iterate step
heredoc
  exit
end

inpfilepath = ARGV[0]
outfilepath = ARGV[1]

albedo = ARGV[2].to_f
ref1   = ARGV[3].to_f
ref2   = ARGV[4].to_f

iterplotfilename = params['iterplot']

if File.exist?(File.expand_path(inpfilepath)) then
  datafile = open(File.expand_path(inpfilepath), 'r+b')
  numrecord = datafile.size / 16
  lut = NArray.sfloat(4, numrecord)
  tmpstr = datafile.read(4*4*numrecord)
  lut[true,true] = NArray.to_na(tmpstr, "sfloat", 4, numrecord)
else
  puts "Error: There's no file."
  exit
end

# call MICO-RE
micore_stdout, micore_status = Open3.capture2("#{MICORE_BIN} #{inpfilepath} #{albedo} #{ref1} #{ref2}")

estimated_tau  = nil
estimated_cder = nil
estimated_cost = nil
niter = nil
iter_taus = []
iter_cders = []
iter_costs = []
iter_ref1s = []
iter_ref2s = []

micore_stdout.each_line do |line|
  if line =~ /^\s*TAU:/ then
    estimated_tau = line.gsub(/^\s*TAU:\s*/,'').to_f
  elsif line =~ /^\s*CDER:/ then
    estimated_cder = line.gsub(/^\s*CDER:\s*/,'').to_f
  elsif line =~ /^\s*COST:/ then
    estimated_cost = line.gsub(/^\s*COST:\s*/,'').to_f
  elsif line =~ /^\s*#\s*iterate step:/ then
    niter = line.gsub(/^\s*#\s*iterate step:\s*/,'').to_i
  elsif line =~ /^\s*#\s*estimated TAU:/ then
    iter_taus.push(line.gsub(/^\s*#\s*estimated TAU:\s*/,'').to_f)
  elsif line =~ /^\s*#\s*estimated CDER:/ then
    iter_cders.push(line.gsub(/^\s*#\s*estimated CDER:\s*/,'').to_f)
  elsif line =~ /^\s*#\s*COST:/ then
    iter_costs.push(line.gsub(/^\s*#\s*COST:\s*/,'').to_f)
  elsif line =~ /^\s*#\s*estimated REF1:/ then
    iter_ref1s.push(line.gsub(/^\s*#\s*estimated REF1:\s*/,'').to_f)
  elsif line =~ /^\s*#\s*estimated REF2:/ then
    iter_ref2s.push(line.gsub(/^\s*#\s*estimated REF2:\s*/,'').to_f)
  end
end

lut[2, true] = lut[2, true] + albedo
lut[3, true] = lut[3, true] + albedo

unqtaus  = lut.transpose(1,0).to_a[0].uniq
unqcders = lut.transpose(1,0).to_a[1].uniq

Gnuplot.open do |gp|
  Gnuplot::Plot.new(gp) do |plot|
    case outfilepath
    when /^none$/
      plot.terminal 'x11'
    when /\.tex$/
      plot.terminal "tikz"
      plot.output File.expand_path(outfilepath)
    when /\.pdf$/
      plot.terminal "pdf"
      plot.output File.expand_path(outfilepath)
    when /\.png$/
      plot.terminal "pngcairo"
      plot.output File.expand_path(outfilepath)
    else
      puts "Error: unknown filetype"
      exit
    end

    plot.title "#{inpfilepath.gsub(/^.*\//,'').gsub(/\.bin$/, '').gsub(/_/, ' ')}"
    plot.xlabel "Reflectance 1"
    plot.ylabel "Reflectance 2"

    plot.xrange "[0.0:#{lut.transpose(1,0).to_a[2].max + 0.1}]"
    plot.yrange "[0.0:#{lut.transpose(1,0).to_a[3].max + 0.1}]"

    plot.set "key off"

    ilbl = 1

    lblcders = [4.0, 9.0, 15.0, 20.0, 25.0, 32.0]
    lbltaus  = [1.0, 5.0, 10.0, 50.0, 100.0]

    unqcders.each do |cder|
      ls = lut.to_a.select do |l|
        l[1] == cder
      end

      pltref1 = []
      pltref2 = []

      ls.each do |l|
        pltref1.push(l[2])
        pltref2.push(l[3])
      end

      if lblcders.include?(cder) then
        plot.set "label #{ilbl} at first #{pltref1[-1] + 0.02},#{pltref2[-1]} '#{cder}'"
        ilbl += 1
      end

      rdclr = 14  + (215 / (unqcders.max - unqcders.min) * (cder - unqcders.min))
      grclr = 191 - (158 / (unqcders.max - unqcders.min) * (cder - unqcders.min))

      plot.data << Gnuplot::DataSet.new([pltref1, pltref2]) do |ds|
        ds.with = "lines"
        ds.linecolor = "rgb \"##{format('%02x', rdclr)}#{format('%02x', grclr)}00\""
      end
    end
    unqtaus.each do |tau|
      ls = lut.to_a.select do |l|
        l[0] == tau
      end

      pltref1 = []
      pltref2 = []

      ls.each do |l|
        pltref1.push(l[2])
        pltref2.push(l[3])
      end

      if lbltaus.include?(tau) then
        plot.set "label #{ilbl} at first #{pltref1[0] - 0.02},#{pltref2[0] + 0.03} '#{tau}'"
        ilbl += 1
      end

      plot.data << Gnuplot::DataSet.new([pltref1, pltref2]) do |ds|
        ds.with = "lines"
        if lbltaus.include?(tau) then
          ds.linecolor = "rgb \"#45a1cf\""
        else
          ds.linecolor = "rgb \"#999999\""
        end
      end
    end

    if ref1 and ref2 then
      # plot observed reflectances
      plot.data << Gnuplot::DataSet.new([[ref1], [ref2]]) do |ds|
        ds.with = "points pointtype 1"
        ds.linecolor = "rgb \"#000000\""
      end

      # estimated tau, cder and cost
      plot.set "label #{ilbl + 1} at first 0.01,#{lut.transpose(1,0).to_a[3].max + 0.08}  'TAU:  #{estimated_tau}'"
      plot.set "label #{ilbl + 2} at first 0.01,#{lut.transpose(1,0).to_a[3].max + 0.03} 'CDER: #{estimated_cder}'"
      plot.set "label #{ilbl + 3} at first 0.01,#{lut.transpose(1,0).to_a[3].max - 0.02} 'COST: #{estimated_cost}'"
    end
  end
end


if iterplotfilename

  combs = [
    ['tau', iter_taus],
    ['cder', iter_cders],
    ['cost', iter_costs],
    ['ref1', iter_ref1s],
    ['ref2', iter_ref2s],
  ]

  combs.each do |varname, itervals|
    Gnuplot.open do |gp|
      Gnuplot::Plot.new(gp) do |plot|
        case iterplotfilename
        when /^none$/
          plot.terminal 'x11'
        when /\.tex$/
          plot.terminal "tikz"
          plot.output (File.expand_path(iterplotfilename.gsub(/\.tex$/, ".#{varname}.tex")))
        when /\.pdf$/
          plot.terminal "pdf"
          plot.output (File.expand_path(iterplotfilename.gsub(/\.pdf$/, ".#{varname}.pdf")))
        when /\.png$/
          plot.terminal "pngcairo"
          plot.output (File.expand_path(iterplotfilename.gsub(/\.png$/, ".#{varname}.png")))
        else
          STDERR.puts "Error: unknown filetype"
          exit
        end

        plot.title "#{varname}"
        plot.xlabel "iterate step"
        plot.ylabel "#{varname}"

        plot.xrange "[1:#{niter}]"
        plot.yrange "[0.0:#{itervals.max + 0.1 * itervals.max}]"

        if varname == 'cost'
          plot.yrange "[#{0.1 * itervals.min}:#{10.0 * itervals.max}]"
          plot.set 'logscale y'
        end

        plot.set "key off"
        plot.set "size square"

        plot.data << Gnuplot::DataSet.new([[*1..niter], itervals]) do |ds|
          ds.with = "lines"
        end
      end
    end
  end
end

