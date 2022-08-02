/*
 * Find dots (PML) in outlined Cells 
 * Measure integrated intensity, nb of dots per cells detect if cell infected 
 * Author Philippe Mailly
 */



import ij.*;
import ij.io.FileSaver;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.MeasureFeret;
import mcib3d.geom2.measurements.MeasureSurface;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;


public class Coloc_genes_2DCells implements PlugIn {

    private final boolean canceled = false;
    private String imageDir = "";
    public static String outDirResults = "";
    private BufferedWriter cellsNumberResults;
    private BufferedWriter cellsAreaResults;
    
     private Genes_2DCells_Processing genes = new Genes_2DCells_Processing(); 

    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        try {
            imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }
            // Find images files
            String fileExt = genes.findImageType(new File(imageDir));
            ArrayList<String> imageFiles = genes.findImages(imageDir, fileExt);
            if (imageFiles == null) {
                IJ.showMessage("Error", "No images found with "+fileExt+" extension");
                return;
            }
            // create output folder
            outDirResults = imageDir + File.separator+ "Results"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            
            // Write headers results for results file
            FileWriter fileResultsNumber = new FileWriter(outDirResults + "CellsNumberResults.xls", false);
            cellsNumberResults = new BufferedWriter(fileResultsNumber);
            try {
                cellsNumberResults.write("ImageName\t#Dapi\t#Gene1\t#Gene2\t#Gene3\t#Dapi/Gene1\t#Dapi/Gene2\t#Dapi/Gene3\t#Gene1/Gene2\t#Gene1/Gene3"
                        + "\t#Gene2/Gene3\t#Gene1/Gene2/Gene3\n");
                cellsNumberResults.flush();
            }
            catch (IOException ex) {
                Logger.getLogger(Coloc_genes_2DCells.class.getName()).log(Level.SEVERE, null, ex);
            }
            FileWriter fileAreaResults = new FileWriter(outDirResults + "CellsAreaResults.xls");
            cellsAreaResults = new BufferedWriter(fileAreaResults);
            try {
                cellsAreaResults.write("ImageName\tDapi area\tGene1 area\tGene2 area\tGene3 area\tDapi/Gene1(area)\tDapi/Gene2(area)\tDapi/Gene3(area)\t"
                        + "Gene1(area)/Gene2\tGene1(area)/Gene3\tGene2(area)/Gene3\tGen1(area)/Gene2)/Gene3\n");
                cellsAreaResults.flush();    
            } catch (IOException ex) {
                Logger.getLogger(Coloc_genes_2DCells.class.getName()).log(Level.SEVERE, null, ex);
            }
            
            
             // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            // Find calibration
            reader.setId(imageFiles.get(0));
            genes.cal = genes.findImageCalib(meta);
            
            // Channels dialog
            
            boolean cancel = genes.dialog();
            if (cancel) {
                IJ.showStatus("Plugin cancelled");
                return;
            }
            // Read images
            for (String f : imageFiles) {
                reader.setId(f);
                String rootName = FilenameUtils.getBaseName(f);
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                options.setQuiet(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                
                // find dapi and genes channel images

                // Open DAPI Channel and detect nucleus
                ImagePlus imgDapi = BF.openImagePlus(options)[0];
                Objects3DIntPopulation dapiPop = genes.stardistCellsPop(imgDapi, 0);
                System.out.println("Total nucleus found = "+dapiPop.getNbObjects());

                // Open gene1 image
                ImagePlus imgGene1 = BF.openImagePlus(options)[1];
                Objects3DIntPopulation gene1Pop = genes.stardistCellsPop(imgGene1, 1);
                System.out.println("Total gene1 = "+gene1Pop.getNbObjects());
                Objects3DIntPopulation gene1FilterPop = genes.geneIntensityFilter(gene1Pop,imgGene1,genes.gene1IntTh);
                System.out.println("Total gene1 after intensity filter = "+gene1FilterPop.getNbObjects());
                genes.closeImages(imgGene1);
                gene1Pop = null;
                
                // Open Gene2 image
                ImagePlus imgGene2 = BF.openImagePlus(options)[2];
                Objects3DIntPopulation gene2Pop = genes.stardistCellsPop(imgGene2, 2);
                System.out.println("Total gene2 = "+gene2Pop.getNbObjects());
                Objects3DIntPopulation gene2FilterPop = genes.geneIntensityFilter(gene2Pop,imgGene2,genes.gene2IntTh);
                System.out.println("Total gene2 after intensity filter = "+gene2FilterPop.getNbObjects());
                genes.closeImages(imgGene2);
                gene2Pop = null;
                
                // Open Gene3 image
                ImagePlus imgGene3  = BF.openImagePlus(options)[3];
                Objects3DIntPopulation gene3Pop = genes.stardistCellsPop(imgGene3, 3);
                System.out.println("Total gene3 = "+gene3Pop.getNbObjects());
                Objects3DIntPopulation gene3FilterPop = genes.geneIntensityFilter(gene3Pop,imgGene3,genes.gene3IntTh);
                System.out.println("Total gene3 after intensity filter = "+gene3FilterPop.getNbObjects());
                genes.closeImages(imgGene3);
                gene3Pop = null;
                
                // Save objects image
                ImageHandler imhDapi = ImageHandler.wrap(imgDapi).createSameDimensions();
                ImageHandler imhGene1 = imhDapi.createSameDimensions();
                ImageHandler imhGene2 = imhDapi.createSameDimensions();
                ImageHandler imhGene3 = imhDapi.createSameDimensions();
                // draw objects pop
                dapiPop.drawInImage(imhDapi);
                if (gene1FilterPop.getNbObjects() > 0)
                    gene1FilterPop.drawInImage(imhGene1);
                if (gene2FilterPop.getNbObjects() > 0)
                    gene2FilterPop.drawInImage(imhGene2);
                if (gene3FilterPop.getNbObjects() > 0)
                    gene3FilterPop.drawInImage(imhGene3);
                ImagePlus[] imgColors = {imhGene3.getImagePlus(), imhGene1.getImagePlus(), imhDapi.getImagePlus(),null,imhGene2.getImagePlus()};
                ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
                imgObjects.setCalibration(genes.cal);
                FileSaver ImgObjectsFile = new FileSaver(imgObjects);
                ImgObjectsFile.saveAsTiff(outDirResults+rootName+"_Objects.tif");
                genes.closeImages(imgDapi);

                // Find colocalized cells

                // gene1/dapi
                Objects3DIntPopulation gene1DapiPop = genes.findColoc(gene1FilterPop, dapiPop);
                System.out.println(gene1DapiPop.getNbObjects() + " gene1 with dapi found");
                // gene2/dapi
                Objects3DIntPopulation gene2DapiPop = genes.findColoc(gene2FilterPop, dapiPop);
                System.out.println(gene2DapiPop.getNbObjects() + " gene2 with dapi found");
                // gene3/dapi
                Objects3DIntPopulation gene3DapiPop = genes.findColoc(gene3FilterPop, dapiPop);
                System.out.println(gene3DapiPop.getNbObjects() + " gene3 with dapi found");
                // gene1 / gene2
                Objects3DIntPopulation gene1Gene2Pop = genes.findColoc(gene1DapiPop, gene2DapiPop);
                System.out.println(gene1Gene2Pop.getNbObjects() + " gene1 with gene2 found");
                // gene1 / gene3
                Objects3DIntPopulation gene1Gene3Pop = genes.findColoc(gene1DapiPop, gene3DapiPop);
                System.out.println(gene1Gene3Pop.getNbObjects() + " gene1 with gene3 found");
                // gene2 / gene3
                Objects3DIntPopulation gene2Gene3Pop = genes.findColoc(gene2DapiPop, gene3DapiPop);
                System.out.println(gene2Gene3Pop.getNbObjects() + " gene2 with gene3 found");
                // gene1 / gene2 / gene3
                Objects3DIntPopulation gene1Gene2Gene3Pop = genes.findColoc3(gene1DapiPop, gene2DapiPop, gene3DapiPop);
                System.out.println(gene1Gene2Gene3Pop.getNbObjects() + " gene1 with gene2 and gene3 found");

                // write results
                IJ.showStatus("Writing results ...");
                cellsNumberResults.write(rootName+"\t"+dapiPop.getNbObjects()+"\t"+gene1FilterPop.getNbObjects()+"\t"+gene2FilterPop.getNbObjects()+"\t"+gene3FilterPop.getNbObjects()+"\t"+
                        gene1DapiPop.getNbObjects()+"\t"+gene2DapiPop.getNbObjects()+"\t"+gene3DapiPop.getNbObjects()+"\t"+gene1Gene2Pop.getNbObjects()+"\t"+gene1Gene3Pop.getNbObjects()+"\t"+
                        gene2Gene3Pop.getNbObjects()+"\t"+gene1Gene2Gene3Pop.getNbObjects()+"\n");
                cellsNumberResults.flush();
                String dapiVol, gene1Vol, gene2Vol, gene3Vol, gene1DapiVol, gene2DapiVol, gene3DapiVol, gene1Gene2Vol, gene1Gene3Vol, gene2Gene3Vol, gene1Gene2Gene3Vol;
                for(Object3DInt nuc : dapiPop.getObjects3DInt()) {
                    float i = nuc.getLabel();
                    dapiVol = new MeasureVolume(nuc).getVolumeUnit() +"\t";
                    gene1Vol = (gene1FilterPop.getObjectByLabel(i)!= null) ? new MeasureVolume(gene1FilterPop.getObjectByLabel(i)).getVolumeUnit() +"\t" : "NaN\t";
                    gene2Vol = (gene2FilterPop.getObjectByLabel(i)!= null) ? new MeasureVolume(gene2FilterPop.getObjectByLabel(i)).getVolumeUnit() +"\t" : "NaN\t";
                    gene3Vol = (gene3FilterPop.getObjectByLabel(i)!= null) ? new MeasureVolume(gene3FilterPop.getObjectByLabel(i)).getVolumeUnit() +"\t" : "NaN\t";
                    gene1DapiVol = (gene1DapiPop.getObjectByLabel(i)!= null) ? new MeasureVolume(gene1DapiPop.getObjectByLabel(i)).getVolumeUnit() +"\t" : "NaN\t";
                    gene2DapiVol = (gene2DapiPop.getObjectByLabel(i)!= null) ? new MeasureVolume(gene2DapiPop.getObjectByLabel(i)).getVolumeUnit() +"\t" : "NaN\t";
                    gene3DapiVol = (gene3DapiPop.getObjectByLabel(i)!= null) ? new MeasureVolume(gene3DapiPop.getObjectByLabel(i)).getVolumeUnit() +"\t" : "NaN\t";
                    gene1Gene2Vol = (gene1Gene2Pop.getObjectByLabel(i)!= null) ? new MeasureVolume(gene1Gene2Pop.getObjectByLabel(i)).getVolumeUnit() +"\t" : "NaN\t";
                    gene1Gene3Vol = (gene1Gene3Pop.getObjectByLabel(i)!= null) ? new MeasureVolume(gene1Gene3Pop.getObjectByLabel(i)).getVolumeUnit() +"\t" : "NaN\t";
                    gene2Gene3Vol = (gene2Gene3Pop.getObjectByLabel(i)!= null) ? new MeasureVolume(gene2Gene3Pop.getObjectByLabel(i)).getVolumeUnit() +"\t" : "NaN\t";
                    gene1Gene2Gene3Vol = (gene1Gene2Gene3Pop.getObjectByLabel(i)!= null) ? new MeasureVolume(gene1Gene2Gene3Pop.getObjectByLabel(i)).getVolumeUnit() +"\t" : "NaN\t";
                    cellsAreaResults.write(rootName+"\t"+dapiVol+gene1Vol+gene2Vol+gene3Vol+gene1DapiVol+gene2DapiVol+gene3DapiVol+gene1Gene2Vol+gene1Gene3Vol+
                            gene2Gene3Vol+gene1Gene2Gene3Vol+"\n");
                    cellsAreaResults.flush();
                }
                   
            }
            
            
                    
        } catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(Coloc_genes_2DCells.class.getName()).log(Level.SEVERE, null, ex);
        }
        try {
            cellsNumberResults.close();
            cellsAreaResults.close();
        } catch (IOException ex) {
            Logger.getLogger(Coloc_genes_2DCells.class.getName()).log(Level.SEVERE, null, ex);
        }
            IJ.showStatus("Process done");
    }
}
