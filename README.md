# FoilSim3JavaStandalone
Port of various NASA projects from a Java applets to standalone Java applications.

FoilSim III Student Version 1.5a
https://www.grc.nasa.gov/www/k-12/airplane/foil3.html

TunnelSys - Tunnel Test Applet Version 1.0a
https://www.grc.nasa.gov/www/k-12/airplane/tunwtest.html

Interactive Wright 1901 Wind Tunnel
https://wright.nasa.gov/airplane/tunnlint.html

While it's easy to add a main method to an applet, and add the applet panel to a Frame, and thus make it runnable and displayable as an application, this is not generally sufficient. The reason is that browsers provide some extra infrastucture for applets to use. To help with this, Jef Poskanzer has written a very useful adapter class called MainFrame which helps provide this infrastructure. 
