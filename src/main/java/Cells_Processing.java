
import Cellpose.CellposeSegmentImgPlusAdvanced;
import Cellpose.CellposeTaskSettings;
import StardistOrion.StarDist2D;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.ZProjector;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DComputation;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Object3DIntLabelImage;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationColocalisation;
import mcib3d.image3d.ImageByte;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.processing.BinaryMorpho;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import org.apache.commons.io.FilenameUtils;



/**
 *
 * @author phm
 */

public class Cells_Processing {
    
    // min size for nucleus µm3
    private double minNuc = 10;
    // max size for nucleus µm3
    private double maxNuc = 1000;
    
    // max size for cells µm3
    private double maxCC1 = maxNuc*2;
    private double maxGFP = maxCC1;
    
    
    // nucleus dilatation
    public float nucDil = 1;
    
    private double nucProbthr = 0.40;
    private double nucOver = 0.25;
    
    // gene Intensity threshold
    public double gfpIntTh = 500;
    public double cc1IntTh = gfpIntTh;
    
    // StarDist
    private Object syncObject = new Object();
    private final double stardistPercentileBottom = 0.2;
    private final double stardistPercentileTop = 99.8;
    private File modelsPath = new File(IJ.getDirectory("imagej")+File.separator+"models");
    private String stardistModel = "";
    private String stardistOutput = "Label Image"; 
    
    // Cellpose
    public int gfpCellPoseDiameter = 50;
    public int cc1CellPoseDiameter = 50;
    private boolean useGpu = true;
    private String[] cellposeModels = {"cyto","nuclei","tissuenet","livecell", "cyto2", "general","CP", "CPx", "TN1", "TN2", "TN3", "LC1",
        "LC2", "LC3", "LC4"};
    public String gfpCellModel = "";
    public String cc1CellModel = "";
    private String cellPoseEnvDirPath = "/home/phm/.conda/envs/cellpose";
    
    // cell association
    private double minColoc = 0.1;     
    private double maxBB = 0;
    
    // pixel size
    private double pixelSize = 0.258;
    public Calibration cal = new Calibration();
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));

    public Boolean canceled = false;
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
    
    /*
    Find starDist models in Fiji models folder
    */
    public String[] findStardistModels() {
        FilenameFilter filter = (dir, name) -> name.endsWith(".zip");
        File[] modelList = modelsPath.listFiles(filter);
        String[] models = new String[modelList.length];
        for (int i = 0; i < modelList.length; i++) {
            models[i] = modelList[i].getName();
        }
        Arrays.sort(models);
        return(models);
    }
    
    /**
     * Find image calibration
     */
    public Calibration findImageCalib(IMetadata meta) {
        Calibration cal = new Calibration();  
        // read image calibration
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);
        return(cal);
    }
    
    /**
     * Find channels name
     * @param imageName
     * @param imageExt
     */
    public ArrayList<String> findChannels (String imageName) throws DependencyException, ServiceException, FormatException, IOException {
        ArrayList<String> channels = new ArrayList<>();
        // create OME-XML metadata store of the latest schema version
        ServiceFactory factory;
        factory = new ServiceFactory();
        OMEXMLService service = factory.getInstance(OMEXMLService.class);
        IMetadata meta = service.createOMEXMLMetadata();
        ImageProcessorReader reader = new ImageProcessorReader();
        reader.setMetadataStore(meta);
        reader.setId(imageName);
        int chs = reader.getSizeC();
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                String channelsID = meta.getImageName(0);
                for (String c : channelsID.replace("_", "-").split("/")) {
                    channels.add(c);
                };
                break;

            case "lif" :
                String[] ch = new String[chs];
                if (chs > 1) {
                    for (int n = 0; n < chs; n++) 
                        if (meta.getChannelExcitationWavelength(0, n) == null)
                            channels.add(Integer.toString(n));
                        else 
                            channels.add(meta.getChannelExcitationWavelength(0, n).value().toString());
                }
                break;
            default :
                chs = reader.getSizeC();
                for (int n = 0; n < chs; n++)
                    channels.add(Integer.toString(n));
        }
        return(channels);         
    }
    
     /**
     * Do Z projection
     * @param img
     * @param projection parameter
     * @return 
     */
    private ImagePlus doZProjection(ImagePlus img, int param) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(param);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }
    
    /**
    * Find background image intensity
    * Z project min intensity
    * read mean intensity
    * @param img 
    */
    public double[] find_background(ImagePlus img) {
      double[] bg = new double[2];
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      ImageProcessor imp = imgProj.getProcessor();
      bg[0] = imp.getStatistics().mean;
      bg[1] = imp.getStatistics().stdDev;
      System.out.println("Background =  " + bg[0] + " +- " + bg[1]);
      closeImages(imgProj);
      return(bg);
    }
    
    /**
     * Dialog
     */
    public ArrayList dialog(List<String> channels, List<String> channelsName) {
        ArrayList ch = new ArrayList();
        String[] models = findStardistModels();
        if (IJ.isWindows())
            cellPoseEnvDirPath = System.getProperty("user.home")+"\\miniconda3\\envs\\CellPose";
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets(0, 120, 0);
        gd.addImage(icon);
        gd.addMessage("Choose channels", new Font(Font.MONOSPACED , Font.BOLD, 12), Color.black);
        int index = 0;
        for (String chName : channelsName) {
            gd.addChoice(chName + ": ", channels.toArray(new String[0]), channels.get(index));
            index++;
        }
        gd.addMessage("Nucleus parameters", new Font(Font.MONOSPACED , Font.BOLD, 12), Color.black);
        gd.addNumericField("Min nucleus vol. (µm3) :", minNuc);
        gd.addNumericField("Max nucleus vol. (µm3) :", maxNuc); 
        gd.addNumericField("Nucleus dilation (µm):", nucDil);
        gd.addMessage("Stardist parameters", Font.getFont("Monospace"), Color.black);
        if (models.length > 0) {
            gd.addChoice("StarDist model :",models, models[0]);
        }
        else {
            gd.addMessage("No StarDist model found in Fiji !!", Font.getFont("Monospace"), Color.red);
            gd.addFileField("StarDist model :", stardistModel);
        }
        gd.addNumericField("Probability : ", nucProbthr, 2);
        gd.addNumericField("Overlap     : ", nucOver, 2);
        gd.addMessage("Cells parameters", new Font(Font.MONOSPACED , Font.BOLD, 12), Color.black);
        gd.addDirectoryField("Cellpose environment path : ", cellPoseEnvDirPath);
        gd.addMessage("GFP cells parameters", Font.getFont("Monospace"), Color.green);
        gd.addChoice("Cellpose model : ", cellposeModels, cellposeModels[4]);
        gd.addNumericField("Min Cell size (µm3)         : ", gfpCellPoseDiameter, 2);
        gd.addNumericField("Min Intensity Threshold :", gfpIntTh, 4);
        gd.addMessage("CC1 parameters", Font.getFont("Monospace"), Color.red);
        gd.addChoice("Cellpose model : ", cellposeModels, cellposeModels[4]);
        gd.addNumericField("Min Cell size (µm3)         : ", cc1CellPoseDiameter, 2);
        gd.addNumericField("Min Intensity Threshold :", cc1IntTh, 4);
        gd.addMessage("Image calibration", new Font(Font.MONOSPACED , Font.BOLD, 12), Color.black);
        gd.addNumericField("Pixel size : ", pixelSize, 3);
        gd.showDialog();
        if (gd.wasCanceled())
            canceled = true;
        for (int i = 0; i < index; i++)
            ch.add(i, gd.getNextChoice());
        if (models.length > 0) {
            stardistModel = modelsPath+File.separator+gd.getNextChoice();
        }
        else {
            stardistModel = gd.getNextString();
        }
        if (stardistModel.isEmpty()) {
            IJ.error("No model specify !!");
            return(null);
        }
        minNuc = (float)gd.getNextNumber();
        maxNuc = (float)gd.getNextNumber();
        nucDil = (float)gd.getNextNumber();
        nucProbthr = gd.getNextNumber();
        nucOver = gd.getNextNumber();
        cellPoseEnvDirPath = gd.getNextString();
        gfpCellModel = gd.getNextChoice();
        cc1CellPoseDiameter = (int)gd.getNextNumber();
        gfpIntTh = gd.getNextNumber();
        cc1CellModel = gd.getNextChoice();
        cc1CellPoseDiameter = (int)gd.getNextNumber();
        cc1IntTh = gd.getNextNumber();
        pixelSize = gd.getNextNumber();
        cal.pixelWidth = cal.pixelHeight = pixelSize;
        return(ch);
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
    
    public Object3DInt DilateObject(Object3DInt obj, float radXYZ) {
        // special case radii = 0
        if ((radXYZ == 0)) {
            return (obj);
        }
        int radXY = (int)(radXYZ/cal.pixelWidth);
        int radZ = (int)(radXYZ/cal.pixelDepth);
        ImageHandler labelImage = new Object3DIntLabelImage(obj).getCroppedLabelImage(radXY, radXY, radZ, 255, false);
        ImageByte binaryMorpho = mcib3d.image3d.processing.BinaryMorpho.binaryMorpho​(labelImage, BinaryMorpho.MORPHO_DILATE,radXY, radZ, 0);
        Object3DInt objDil = new Object3DInt(binaryMorpho);
        objDil.translate(-radXY, -radXY, -radZ);
        objDil.setVoxelSizeXY(obj.getVoxelSizeXY());
        objDil.setVoxelSizeZ(obj.getVoxelSizeZ());
        objDil.setUnit(obj.getUnit());
        objDil.setLabel(obj.getLabel());
        return objDil;
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
    
   
    /** Look for all nuclei
     * work on images <= 1024
     * split image if size > 1024
     * return nuclei population
     */
    public Objects3DIntPopulation stardistCellsPop(ImagePlus img){
        // resize to be in a stardist-friendly scale
        int width = img.getWidth();
        int height = img.getHeight();
        float factor = 0.5f;
        boolean resized = false;
        ImagePlus imgStar = null;
        if (img.getWidth() > 1024) {
            imgStar = img.resize((int)(width*factor), (int)(height*factor), 1, "none");
            resized = true;
        }
        else
            imgStar = new Duplicator().run(img);
        
        IJ.run(imgStar, "Remove Outliers", "block_radius_x=5 block_radius_y=5 standard_deviations=1 stack");
        File starDistModelFile = new File(stardistModel);
        StarDist2D star = new StarDist2D(syncObject, starDistModelFile);
        double probTh = nucProbthr;
        double over = nucOver;
        double min = minNuc;
        double max = maxNuc;
        star.setParams(stardistPercentileBottom, stardistPercentileTop, probTh, over, stardistOutput);
        star.loadInput(imgStar);
        star.run();
        // label in 3D
        ImagePlus imgLabels = (resized) ? star.associateLabels().resize(width, height, 1, "none") : star.associateLabels();
        imgLabels.setCalibration(cal);
        ImageHandler imh = ImageHandler.wrap(imgLabels);
        Objects3DIntPopulation pop = sizeFilterPop(new Objects3DIntPopulation(imh), min, max);
        closeImages(imgLabels);
        closeImages(imgStar);
        return(pop);
    }
    

    
/**
 * Find cells with cellpose
 * @param img
 * @param type
 * @return 
 */
    public Objects3DIntPopulation cellPoseCellsPop(ImagePlus img, String cellModel, int cellPoseDiameter){
        CellposeTaskSettings settings = new CellposeTaskSettings(cellModel, 1, cellPoseDiameter, cellPoseEnvDirPath);
        settings.setStitchThreshold(0.25); 
        settings.setFlowTh(0.6);
        settings.useGpu(true);
        ImagePlus imgIn = new Duplicator().run(img);
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgIn);
        ImagePlus cellpose_img = cellpose.run(); 
        closeImages(imgIn);
        cellpose_img.setCalibration(cal);
        ImageHandler imh = ImageHandler.wrap(cellpose_img);
        Objects3DIntPopulation pop = sizeFilterPop(new Objects3DIntPopulation(imh), minNuc, maxCC1);
        imh.closeImagePlus();
        closeImages(cellpose_img);
        return(pop);
    }
    
    public Objects3DIntPopulation intensityFilter(Objects3DIntPopulation pop, ImagePlus img, double intTh) {
        Objects3DIntPopulation newPop = new Objects3DIntPopulation();
        if (pop.getNbObjects() > 0) {
            ImageHandler imh = ImageHandler.wrap(img);
            for (Object3DInt obj : pop.getObjects3DInt()) {
                double maxInt = new MeasureIntensity(obj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_MAX);
                if (maxInt > intTh)
                    newPop.addObject(obj);
            }
        }
        return(newPop);
    }
    
   /**
    * Find coloc dapi and cells
    * return first population coloc
    */
    public Objects3DIntPopulation findColoc(Objects3DIntPopulation pop1, Objects3DIntPopulation pop2, String ch, ArrayList<Nucleus> nuclei) {
        Objects3DIntPopulation colocPop = new Objects3DIntPopulation();
        if (pop1.getNbObjects() > 0 && pop2.getNbObjects() > 0) {
            MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(pop1, pop2);
            for (Object3DInt obj1 : pop1.getObjects3DInt()) {
                for (Object3DInt obj2 : pop2.getObjects3DInt()) {
                    double colocVal = coloc.getValueObjectsPair(obj1, obj2);
                    if (colocVal > 0) {
                        colocPop.addObject(obj1);
                        Optional<Nucleus> nuc = nuclei.stream().filter(i -> i.getLabel() == obj1.getLabel()).findFirst();
                        switch (ch) {
                                case "GFP" :
                                    nuc.get().setisGFP(true);
                                    break;
                                case "CC1" :
                                    nuc.get().setisCC1(true);
                                    break;    
                        }
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
    * Measure mean intensity of dilated nucleus in NG2 image
    * 
    */
    public Objects3DIntPopulation measureIntensity(Objects3DIntPopulation dapiPop, ImagePlus imgNG2, ArrayList<Nucleus> nuclei, double bg) {
        ImageHandler imh = ImageHandler.wrap(imgNG2);
        double meanInt = 0;
        double meanIntCor = 0;
        float radXY = (float)(nucDil/cal.pixelWidth);
        float radZ = (float)(nucDil/cal.pixelDepth);
        Objects3DIntPopulation  nucDilated = new Objects3DIntPopulation();
        for (Object3DInt obj : dapiPop.getObjects3DInt()) {
            Object3DInt objDil = new Object3DComputation(obj).getObjectDilated(radXY, radXY, radZ);
            nucDilated.addObject(objDil);
        }
        for(Object3DInt obj : nucDilated.getObjects3DInt()) {
            MeasureIntensity intensity = new MeasureIntensity(obj, imh);
            Optional<Nucleus> nuc = nuclei.stream().filter(i -> i.getLabel() == obj.getLabel()).findFirst();
            if (imgNG2 != null) {
                meanInt = intensity.getValueMeasurement(MeasureIntensity.INTENSITY_AVG);
                meanIntCor = meanInt - bg;
            }
            nuc.get().setNG2MeanInt(meanInt);
            nuc.get().setNG2MeanIntCor(meanIntCor);
        }
        return(nucDilated);
    }
    
    
}
    