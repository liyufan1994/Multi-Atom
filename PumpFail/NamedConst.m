classdef NamedConst
   properties (Constant)
      
      % These are pump model constants
      sarr=[5 1 5 14 3 19 1 1 4 22]
      tarr=[94.32 15.72 62.88 125.76 5.24 31.44 1.048 1.048 2.096 10.480]

      alpha=1.82
      gam=0.01
      sig=1
      
      cent0=[0.07 0.15 0.10 0.12 0.62 0.61 0.83 0.82 1.29 1.84 2.47]
      Radi=1*[0.027 0.0921 0.0401 0.0309 0.2934 0.1343 0.5295 0.5319 0.5752 0.3904 0.7189]

      M=10000;
      Mout=10000; 
      dig=6;
   end
end

