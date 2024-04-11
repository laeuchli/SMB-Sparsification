function setRandomSeed(s)
%SETRANDOMSEED Fix a random seed.
%   SETRANDOMSEED(S) sets the random seed of MATLAB to the value S. This
%   helps obtain deterministic results for random-dependent algorithms.
%
%   See also: RAND, RANDSTREAM.

stream = RandStream('mt19937ar','Seed',s);%RandStream('mlfg6331_64','RandnAlg','Polar','Seed', s);
%RandStream.setDefaultStream(stream);

end
