
//import FoilSim.Foil;
//import Tunl.Foil;
import WTest.Wys;
import Acme.MainFrame;

public class Driver {
    public static void main(String[] args){
        FoilSim.Foil foilSimFoil = new FoilSim.Foil();
        Tunl.Foil tunlFoil = new Tunl.Foil();
        WTest.Wys wTestFoil = new WTest.Wys();
        
        MainFrame foilSimFoilMainFrame = new MainFrame(foilSimFoil, 750, 500);
        MainFrame tunlFoilMainFrame = new MainFrame(tunlFoil, 750, 500);
        MainFrame wTestFoilMainFrame = new MainFrame(wTestFoil, 740, 690);
    }
}