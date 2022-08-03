
import Orion_Stardist.StarDist2D;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import javax.swing.ImageIcon;
import loci.formats.meta.IMetadata;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationColocalisation;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import org.apache.commons.io.FilenameUtils;



/**
 *
 * @author phm
 */

public class Genes_2DCells_Processing {
    
    // min size for nucleus µm2
    private double minNuc = 30;
    // max size for nucleus µm2
    private double maxNuc = 300;

    // Stardist parameters
    private final double stardistPercentileBottom = 0.2;
    private final double stardistPercentileTop = 99.8;    
    private String stardistOutput = "Label Image"; 

    private double nucProbthr = 80;
    private double gene1Probthr = 80;
    private double gene2Probthr = 80;
    private double gene3Probthr = 80;
    
    private double nucOver = 0.25;
    private double gene1Over = 0.25;
    private double gene2Over = 0.25;
    private double gene3Over = 0.25;
    
    // gene Intensity threshold
    public double gene1IntTh = 500;
    public double gene2IntTh = gene1IntTh;
    public double gene3IntTh = gene1IntTh;
    
    // pixel size
    private double pixelSize = 0.258;
    public Calibration cal = new Calibration();
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));

    private CLIJ2 clij2 = CLIJ2.getInstance();
    
     /**
     * check  installed modules
     * @return 
     */
    public boolean checkInstalledModules() {
        // check install
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("uk.ac.sussex.gdsc.utils.DifferenceOfGaussians_PlugIn");
        } catch (ClassNotFoundException e) {
            IJ.log("GDSC Suite not installed, please install from update site");
            return false;
        }
        
        return true;
    }
    
    /**
     * Dialog
     */
    public boolean dialog() {
        Boolean canceled = false;
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets(0, 120, 0);
        gd.addImage(icon);
        gd.addMessage("Nucleus parameters", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Probability : ", nucProbthr, 2);
        gd.addNumericField("Overlap     : ", nucOver, 2);
        gd.addMessage("Genes parameters", Font.getFont("Monospace"), Color.DARK_GRAY);
        gd.addMessage("488 Genes", Font.getFont("Monospace"), Color.green);
        gd.addNumericField("Probability : ", gene1Probthr, 2);
        gd.addNumericField("Overlap     : ", gene1Over, 2);
        gd.addNumericField("Max Intensity Threshold :", gene1IntTh, 4);
        gd.addMessage("561 Genes", Font.getFont("Monospace"), Color.red);
        gd.addNumericField("Probability : ", gene2Probthr, 2);
        gd.addNumericField("Overlap     : ", gene2Over, 2);
        gd.addNumericField("Max Intensity Threshold :", gene2IntTh, 4);
        gd.addMessage("642 Genes", Font.getFont("Monospace"), Color.MAGENTA);
        gd.addNumericField("Probability : ", gene3Probthr, 2);
        gd.addNumericField("Overlap     : ", gene3Over, 2);
        gd.addNumericField("Max Intensity Threshold :", gene3IntTh, 4);
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.DARK_GRAY);
        gd.addNumericField("Pixel size : ", pixelSize, 3);
        gd.showDialog();
        if (gd.wasCanceled())
            canceled = true;
        String imgFolder = gd.getNextString();
        nucProbthr = gd.getNextNumber();
        nucOver = gd.getNextNumber();
        gene1Probthr = gd.getNextNumber();
        gene1Over = gd.getNextNumber();
        gene1IntTh = gd.getNextNumber();
        gene2Probthr = gd.getNextNumber();
        gene2Over = gd.getNextNumber();
        gene2IntTh = gd.getNextNumber();
        gene3Probthr = gd.getNextNumber();
        gene3Over = gd.getNextNumber();
        gene3IntTh = gd.getNextNumber();
        pixelSize = gd.getNextNumber();
        cal.pixelWidth = cal.pixelHeight = pixelSize;
        return(canceled);
    }
    
     /**
     * Find image type
     */
    public String findImageType(File imagesFolder) {
        String ext = "";
        String[] files = imagesFolder.list();
        for (String name : files) {
            String fileExt = FilenameUtils.getExtension(name);
            switch (fileExt) {
               case "nd" :
                   ext = fileExt;
                   break;
                case "czi" :
                   ext = fileExt;
                   break;
                case "lif"  :
                    ext = fileExt;
                    break;
                case "isc2" :
                    ext = fileExt;
                    break;
            }
        }
        return(ext);
    }
    
    
   /**
     * Find images in folder
     */
    public ArrayList findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No Image found in "+imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equalsIgnoreCase(imageExt))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
    
     /**
     * Find image calibration
     * @param meta
     * @return 
     */
    public Calibration findImageCalib(IMetadata meta) {
        cal = new Calibration();  
        // read image calibration
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        return(cal);
    }
    
    
     /**
     * Filter population by size
     */
    public Objects3DIntPopulation sizeFilterPop(Objects3DIntPopulation pop, double volMin, double volMax) {
        Objects3DIntPopulation popF = new Objects3DIntPopulation();
        for (Object3DInt object : pop.getObjects3DInt()) {
            double vol = new MeasureVolume(object).getVolumeUnit();
            if ((vol >= volMin) && (vol <= volMax)) {
                popF.addObject(object);
            }
        }
        popF.setVoxelSizeXY(cal.pixelWidth);
        popF.setVoxelSizeZ(cal.pixelDepth);
        return(popF);
    }
        
    /**
     *
     * @param img
     */
    public void closeImages(ImagePlus img) {
        img.flush();
        img.close();
    }

     /* Median filter 
     * Using CLIJ2
     * @param ClearCLBuffer
     * @param sizeXY
     * @param sizeZ
     */ 
    public ClearCLBuffer median_filter(ClearCLBuffer  imgCL, double sizeXY, double sizeZ) {
        ClearCLBuffer imgCLMed = clij2.create(imgCL);
        clij2.mean3DBox(imgCL, imgCLMed, sizeXY, sizeXY, sizeZ);
        clij2.release(imgCL);
        return(imgCLMed);
    }
    
    
   /**
     * return objects population in an binary image
     * @param img
     * @return pop objects population
     */

    public Objects3DPopulation getPopFromImage(ImagePlus img) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        return pop;
    } 
    
    
    /** Look for all nuclei
     * work on images <= 2048
     * split image if size > 2048
     * return nuclei population
     */
    public Objects3DIntPopulation stardistCellsPop(ImagePlus img, int type){
        // resize to be in a stardist-friendly scale
        int width = img.getWidth();
        int height = img.getHeight();
        float factor = 0.25f;
        boolean resized = false;
        ImagePlus imgStar = null;
        if (img.getWidth() > 1024) {
            imgStar = img.resize((int)(width*factor), (int)(height*factor), 1, "none");
            resized = true;
        }
        else
            imgStar = new Duplicator().run(img);
        
        IJ.run(imgStar, "Remove Outliers", "block_radius_x=10 block_radius_y=10 standard_deviations=1 stack");
        StarDist2D star = new StarDist2D();
        double probTh = 0;
        double over = 0;
        double min = 0;
        double max = 0;
        switch (type) {
            case 0 :
                probTh = nucProbthr;
                over = nucOver;
                min = minNuc;
                max = maxNuc;
                break;
            case 1 : 
                probTh = gene1Probthr;
                over = gene1Over;
                min = 0;
                max = Double.MAX_VALUE;
                break; 
            case 2 : 
                probTh = gene2Probthr;
                over = gene2Over;
                min = 0;
                max = Double.MAX_VALUE;
                break; 
            case 3 : 
                probTh = gene3Probthr;
                over = gene3Over;
                min = 0;
                max = Double.MAX_VALUE;
                break; 
        } 
        
        star.setParams(stardistPercentileBottom, stardistPercentileTop, probTh, over, stardistOutput);
        star.checkImgSize(imgStar);
        star.loadInput(imgStar);
        star.run();
        Img<? extends RealType<?>> img1 = star.label.getImgPlus().getImg();
        ImagePlus imgLab = (resized) ? ImageJFunctions.wrap((RandomAccessibleInterval)img1, "Labelled").resize(width, height, 1, "none") : ImageJFunctions.wrap((RandomAccessibleInterval)img1, "Labelled");
        
        imgLab.setCalibration(cal);
        Objects3DIntPopulation nucPop = new Objects3DIntPopulation(ImageHandler.wrap(imgLab));
        Objects3DIntPopulation popFilter = sizeFilterPop(nucPop, min, max);
        closeImages(imgLab);
        closeImages(imgStar);
        return(popFilter);
    }

    
    
    public Objects3DIntPopulation geneIntensityFilter(Objects3DIntPopulation pop, ImagePlus img, double intTh) {
        Objects3DIntPopulation newPop = new Objects3DIntPopulation();
        ImageHandler imh = ImageHandler.wrap(img);
        for (Object3DInt obj : pop.getObjects3DInt()) {
            double maxInt = new MeasureIntensity(obj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_MAX);
            if (maxInt > intTh)
                newPop.addObject(obj);
        }
        return(newPop);
    }
    
   /**
    * Find coloc with 2 genes
    * return first population coloc
    */
    public Objects3DIntPopulation findColoc(Objects3DIntPopulation pop1, Objects3DIntPopulation pop2) {
        Objects3DIntPopulation colocPop = new Objects3DIntPopulation();
        if (pop1.getNbObjects() > 0 && pop2.getNbObjects() > 0) {
            MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(pop1, pop2);
            //coloc.getResultsTableOnlyColoc().show("coloc");
            for (Object3DInt obj1 : pop1.getObjects3DInt()) {
                for (Object3DInt obj2 : pop2.getObjects3DInt()) {
                    double colocVal = coloc.getValueObjectsPair(obj1, obj2);
                    if (colocVal > 0) {
                        colocPop.addObject(obj1); 
                        break;
                    }
                }
            }
            colocPop.setVoxelSizeXY(cal.pixelWidth);
            colocPop.setVoxelSizeZ(cal.pixelDepth);
        }
        return(colocPop);
    }
    
/**
    * Find coloc with 3 genes
    * return first population coloc
    */
    public Objects3DIntPopulation findColoc3(Objects3DIntPopulation pop1, Objects3DIntPopulation pop2, Objects3DIntPopulation pop3) {
        Objects3DIntPopulation   coloc1 = findColoc(pop1, pop2);
        Objects3DIntPopulation   coloc2 = findColoc(coloc1, pop3);
       return(coloc2);
    }
}
    