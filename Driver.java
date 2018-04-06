
import Acme.MainFrame;

public class Driver {
    public static void main(String[] args){
        FoilSim.Foil foilSim = new FoilSim.Foil();
        Tunl.Foil tunl = new Tunl.Foil();
        WTest.Wys wTest = new WTest.Wys();
        Mach.Mach mach = new Mach.Mach();
        Isentrop.Isentrop isentrop = new Isentrop.Isentrop();
        Shock.Shock shock = new Shock.Shock();
        Shockc.Shock shockc = new Shockc.Shock();
        Mshock.Mshock mshock = new Mshock.Mshock();
        Moc.Moc moc = new Moc.Moc();
        Nozzle.Nozzle nozzle = new Nozzle.Nozzle();
        Sfs.Moc sfs = new Sfs.Moc();
        EngineSimU.Turbo turbo = new EngineSimU.Turbo();
        FoilSimU.Foil foilSimU = new FoilSimU.Foil();
        
        MainFrame foilSimMainFrame = new MainFrame(foilSim, 750, 500);
        MainFrame tunlMainFrame = new MainFrame(tunl, 750, 500);
        MainFrame wTestMainFrame = new MainFrame(wTest, 740, 690);
        MainFrame machFoilMainFrame = new MainFrame(mach, 440, 390);
        MainFrame isentropMainFrame = new MainFrame(isentrop, 630, 345);
        MainFrame shockMainFrame = new MainFrame(shock, 645, 435);
        MainFrame shockcMainFrame = new MainFrame(shockc, 750, 500);
        MainFrame mshockMainFrame = new MainFrame(mshock, 750, 500);
        MainFrame mocMainFrame = new MainFrame(moc, 750, 820);
        MainFrame nozzleMainFrame = new MainFrame(nozzle, 750, 500);
        MainFrame sfsMainFrame = new MainFrame(sfs, 750, 550);
        MainFrame turboMainFrame = new MainFrame(turbo, 750, 550);
        MainFrame foilSimUMainFrame = new MainFrame(foilSimU, 750, 550);
    }
}