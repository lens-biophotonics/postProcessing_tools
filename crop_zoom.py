'''
Libraries
'''

import tifffile as tiff
from scipy.ndimage import zoom
from tqdm import tqdm
import concurrent.futures
import os
import argparse
import yaml

'''
Helper Functions
'''


def get_omeTifs_list(directory_path):
    """
    Retrieves a list of all files in the specified directory with the `.ome.tif` extension.

    Args:
        directory_path (str): Path to the directory containing OME-TIF files.

    Returns:
        list: A list of file names that have the `.ome.tif` extension.
    """
    # List containing .ome.tif file names in the given directory.
    omeTif_list = [f for f in os.listdir(directory_path) if f.endswith(".ome.tif")]

    return omeTif_list


def read_omeTif(path):
    """
    Reads an OME-TIF file and extracts its metadata and image data.

    Args:
        path (str): Path to the OME-TIF file.

    Returns:
        tuple:
            - omeTif_metadata (str): OME metadata as an XML string.
            - omeTif_imgData (numpy.ndarray): Image data from the OME-TIF file as a NumPy array.
    """
    # Open the OME-TIF file.
    with tiff.TiffFile(path) as tif:
        omeTif_metadata = tif.ome_metadata
        # Load the image data and convert it to a NumPy array
        omeTif_imgData = tif.asarray()

    return (omeTif_metadata, omeTif_imgData)


def save_transformed_omeTif(output_filePath, omeTif_imgData, omeTif_metadata):
    """
    Saves OME-TIFF file at the given path.

    Args:
        output_filePath (str): Path to sav the OME-TIF file.

        omeTif_imgData (numpy.ndarray):  Image data of the OME-TIF file.

        omeTif_metadata (str): OME-TIF file metadata
    """
    # Open the TIFF writer in BigTIFF format to handle large files
    with tiff.TiffWriter(output_filePath, bigtiff=True) as tif:
        # Write the OME-TIF file
        tif.write(
            omeTif_imgData,  # OME-TIF image data
            photometric="minisblack",  # Grayscale image specification
            metadata={"axes": "ZYX", "ome": omeTif_metadata},  # Metadata includes axes and OME information
            compression="zlib",  # Use zlib compression to reduce file size
        )


def crop_zoom_ometif(file_name, directory_path, output_dirPath, x_axis):
    """
    Crops the OME-TIF file to the desired x dimensions and then 
    zooms it back to the original X shape.

    Args:
        file_name (str): Name of the OME-TIFF file.

        directory_path (str): Path to the directory containing the OME-TIF file.

        output_dirPath (str): Path to the directory to sav the OME-TIF file.

        x_axis (tuple): A tuple with X1 and X2 coordinates for cropping.
    """
    # creats the path to/save the OME-TIF file.
    file_path = os.path.join(directory_path, file_name)
    output_filePath = os.path.join(output_dirPath, file_name)
    
    # Read
    omeTif_metadata, omeTif_imgData = read_omeTif(file_path)
    
    # Get the X shape
    _, _, X_shape = omeTif_imgData.shape

    # Crop the image to the shape -> (Z, Y, x2 - x1)
    omeTif_imgData = omeTif_imgData[:, :, x_axis[0]:x_axis[1]] 

    # Zoom factor for each dimension to go back to original shape -> (Z, Y, X)
    zoom_factors = (
        1.0,  # no change to Z
        1.0,  # no change to Y
        X_shape / (x_axis[1] - x_axis[0])  # stretch X dimension back
    )
    # Current shape is (Z, Y, x2 - x1). Original is (Z, Y, X).
    # Zoom with order=1 -> bilinear
    omeTif_imgData = zoom(omeTif_imgData, zoom=zoom_factors, order=1)
    
    # Save
    save_transformed_omeTif(output_filePath, omeTif_imgData, omeTif_metadata)

    return 1


def parallel_crop_zoom(directory_path, output_dirPath, omeTif_list, x_axis, tiles, max_workers=2):
    """
    In Parallel, crop the OME-TIFF files and zoom back to the original dimensions.

    Args:

        directory_path (str): Path to the directory containing the OME-TIF files.

        output_dirPath (str): Path to the directory to sav the OME-TIF file.

        omeTif_list (list): List of `.ome.tif` file names at the directory_path.

        x_axis (tuple): A tuple with X1 and X2 coordinates for cropping.

        tiles (tuple): A tuple to specify the start and end tile to perform the crop and zoom.

        max_worker (int): Number of CPUs working in parallel.
    """
    # counter
    done = 0
    
    # ProcessPoolExecutor for parallel processing 
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all crop_zoom jobs
        futures = [
            executor.submit(
                crop_zoom_ometif, file_name, directory_path, output_dirPath, x_axis
            )
            for file_name in omeTif_list
        ]

        # Collect results with a progress bar
        for future in tqdm(
            concurrent.futures.as_completed(futures), 
            total=len(futures), 
            desc=f"Applying crop and zoom on the set ({tiles[0]}, {tiles[1]}) parallel"
        ):
            # error check
            try:
                done += future.result()
            except Exception as e:
                print(f"Error processing a file: {e}")

    print(f"All {done} files have been processed and saved.")



'''
Main Function
'''


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Crop and zoom OME-TIFF files using parameters from a YAML file.")
    parser.add_argument("-p", type=str, required=True, help="Path to the YAML configuration file.")
    args = parser.parse_args()

    # Load parameters from the YAML file
    with open(args.p, 'r') as file:
        params = yaml.safe_load(file)

    # Directory containing and to save the OME-TIF files
    directory_path = params['directory_path']
    output_dirPath = params['output_dirPath']
    # A list containing a lists of X coordinates under the index [0] and [1] for 
    # cropping and,
    # the starting and ending tiles to perform that cropping on as under the index
    # [2] and [3].
    crop_dict = params['crop_dict']
    # number of cpus to work in parallel
    max_workers = params.get('max_workers', 1)  # Default to 1 if not specified

    # Create output directory if it doesn't exist
    try:
        os.makedirs(output_dirPath, exist_ok=True)
        print(f"Directory {output_dirPath} created successfully")
    except OSError as error:
        print(f"Directory {output_dirPath} cannot be created")


    # Get the list of OME-TIF files
    omeTif_list = get_omeTifs_list(directory_path)
    omeTif_list.sort()

    # for every X coordinates crop
    for group in crop_dict:
        x_axis = (group[0], group[1])
        tiles = (group[2], group[3])
        # get the tiles name to perform the cropping on.
        omeTif_list_temp = [
            omeTif_list[x + tiles[0]]           
            for x in range(tiles[1] - tiles[0] + 1)
        ]
        # perform crop and zoom in parallel
        parallel_crop_zoom(directory_path, output_dirPath, omeTif_list_temp, x_axis, tiles, max_workers)



if __name__ == "__main__":
    main()