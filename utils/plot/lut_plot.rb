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
#

require 'narray'
require 'gnuplot'


if ARGV.size < 3 then
  STDERR.puts <<heredoc
Usage: ruby #{$0} input output numrecord
  input:  lut filename
  output: output png filename
  numrecord = (# of tau) * (# of cder)
heredoc
  exit
end

inpfilepath = ARGV[0]
outfilepath = ARGV[1]
numrecord   = ARGV[2].to_i

if File.exist?(File.expand_path(inpfilepath)) then
  datafile = open(File.expand_path(inpfilepath), 'r+b')
  lut = NArray.sfloat(4, numrecord)
  tmpstr = datafile.read(4*4*numrecord)
  lut[true,true] = NArray.to_na(tmpstr, "sfloat", 4, numrecord)
else
  puts "There's no file."
  exit
end

unqtaus  = lut.transpose(1,0).to_a[0].uniq
unqcders = lut.transpose(1,0).to_a[1].uniq

Gnuplot.open do |gp|
  Gnuplot::Plot.new(gp) do |plot|
    if outfilepath =~ /none/ then
      plot.terminal 'x11'
    elsif outfilepath =~ /\.tex$/ then
      plot.terminal "tikz"
      plot.output File.expand_path(outfilepath)
    elsif outfilepath =~ /\.pdf$/ then
      plot.terminal "pdf"
      plot.output File.expand_path(outfilepath)
    elsif outfilepath =~ /\.png$/ then
      plot.terminal "png"
      plot.output File.expand_path(outfilepath)
    end

    plot.title "#{inpfilepath.gsub(/^.*\//,'').gsub(/\.bin$/, '').gsub(/_/, ' ')}"
    plot.xlabel "Reflectance 1"
    plot.ylabel "Reflectance 2"

    plot.xrange "[0.0:#{lut.transpose(1,0).to_a[2].max + 0.1}]"
    plot.yrange "[0.0:#{lut.transpose(1,0).to_a[3].max + 0.1}]"

    plot.set "key off"

    clr = 0
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

      rdclr = 14  + (215 / unqcders.length * clr)
      grclr = 191 - (158 / unqcders.length * clr)

      plot.data << Gnuplot::DataSet.new([pltref1, pltref2]) do |ds|
        ds.with = "lines"
        ds.linecolor = "rgb \"##{format('%02x', rdclr)}#{format('%02x', grclr)}00\""
      end
      clr += 1
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
        ds.linecolor = "rgb \"#999999\""
      end
    end

  end
end


