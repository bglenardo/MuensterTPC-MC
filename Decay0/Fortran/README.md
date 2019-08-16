DECAY0
========

Changelog:
----------



14/01/2013 :
     - implemented Xe134 2nubb decay (adjusted Q value)
            * Qbb= 0.843 MeV
     - implemented Xe124 ECEC, ECb+, b+b+ (adjusted Q value and K X-ray) 
            * Q = 3.068 MeV
            * K_\alpha X-ray : 27.8 keV
	- Replaced old continuation character (+) by the one for tab form (any digit except 0)
		* sed "/^[     +]*/ s/     +/     9/g" decay0.1.for > decay0.1.new.for

11/01/2013 : 
	- Imported original DECAY0 code and started converting into DECAY0.1 so that Eclipse Photran can parse it successfully.
	- Replaced the old comments (starting with 'c') with the following command:
		* sed "/^[c]*/ s/c/\!c/g" test.txt 
11/06/2019 :
	- Added excited states of Te-124 (ll. 1313ff.)
		!c Added excited nuclear states relevant for double positron decays
		!c Spin and parity of states taken from:
		!c https://doi.org/10.1016/j.nds.2008.06.001
		!c http://dx.doi.org/10.1155/2013/505874
		!c Keep in mind that it is not clear if 2790 keV is 0+, 2+ or 4+
		if(ilevel.eq.0) levelE=0
	        if(ilevel.eq.1) levelE=602.73
	        if(ilevel.eq.2) levelE=1325.51
	        if(ilevel.eq.3) levelE=1675.28
	        if(ilevel.eq.4) levelE=2790.41
	        if(ilevel.eq.0.or.ilevel.eq.3.or.ilevel.eq.4) itrans02=0
	        if(ilevel.eq.1.or.ilevel.eq.2) itrans02=2
	- Added excited Te-124 states to interface (ll. 314ff.) and fixed bug that the interface would list Ba-136 states for Xe-124 decay
		elseif((chn(1:2).eq.'Xe'.or.chn(1:2).eq.'XE'.or.
     9             	chn(1:2).eq.'xe').and.chn(3:5).eq.'124') then
	124	      print *,'124-Te level:   0. 0+ (gs)     {0 MeV}'
		      print *,'                1. 2+ (1)  {0.60273 MeV}'
		      print *,'                2. 2+ (2)  {1.32551 MeV}'
		      print *,'                3. 0+ (1)  {1.65728 MeV}'
		      print *,'                4. 0+ (2)  {2.79041 MeV}'
		      print 1
		      read *,ilevel
		      if(ilevel.lt.0.or.ilevel.gt.4) then
			 print *,'   '
			 go to 124
		      endif
16/08/2019 :
	- Added a nuclear deexcitation subroutine for simulation of Xe-124 ECEC/ECb+/b+b+ decays into excited states of Te-124