import argparse
import prad.screen
import prad.deflect

def get_input_data():

    parser = argparse.ArgumentParser(
                description='This script is used to reconstruct the path-'
                            'integrated magnetic field from a flux image.')

    parser.add_argument('flux_image',
                        action='store', type=str,
                        help='Path to flux image.')

    parser.add_argument('flux_mask',
                        action='store', type=str,
                        help='Path to flux mask.')
    
    parser.add_argument('source_to_plasma',
                        action='store', type=float,
                        help='Distance (cm) from the proton source to the '
                             'plasma region.')
    
    parser.add_argument('plasma_length',
                        action='store', type=float,
                        help='Distance (cm) across the '
                             'plasma region.')
    
    parser.add_argument('plasma_to_screen',
                        action='store', type=float,
                        help='Distance (cm) from the plasma region to the '
                             'detector.')
    
    parser.add_argument('pixel_size',
                        action='store', type=float,
                        help='Pixel size (cm) of the detector bins.')
    
    parser.add_argument('proton_velocity',
                        action='store', type=float,
                        help='Proton velocity (cm/s).')

    parser.add_argument('chk_interval',
                        action='store', type=int,
                        help='Step interval to save checkpoint files.')

    parser.add_argument('save_dir',
                        action='store', type=str,
                        help='Directory to save checkpoint files.')
    
    args = parser.parse_args()

    return args

def reconstruct_Bpath():
    args = get_input_data()
    
    flux_image = np.loadtxt(args.flux_image, delimiter=',')
    flux_mask = np.loadtxt(args.flux_mask, delimiter=',')

    plasma_x, plasma_y, flux0, flux = prad.screen.setup(
                                            args.source_to_plasma,
                                            args.plasma_length,
                                            args.plasma_to_screen,
                                            args.pixel_size,
                                            flux_image,
                                            flux_mask)
    
    # TODO Add functionality to let the user choose their tolerance level &
    # timestep.
    
    phin, phix, phiy = prad.screen.solve(
                            plasma_x, plasma_y,
                            flu0, flux,
                            1e-8, 1,
                            chk=True,
                            interval=args.chk_interval,
                            nan_exception=True,
                            save_dir=args.save_dir)

if __name__=='__main__':
    reconstruct_Bpath()
