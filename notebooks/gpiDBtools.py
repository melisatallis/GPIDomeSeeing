"""
keep helpful wrappers for common GPI Database actions here.
"""
#!/usr/bin/env python
import glob
import os
import sys
import numpy as np
import time
from astropy.io import fits
from datetime import datetime
#import gpifilesdb
#import shutil
#from DB_delete_date import test_entries, delete_entries



def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [[Y]/n] "
    elif default == "no":
        prompt = " [y/[N]] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")



##############################################################################
##############################################################################



def tstamp_from_whenstr(whenstr):
    if len(whenstr) != 14:
        sys.exit('whenstr must be YYYYMMDDhhmmss')

    Y = whenstr[:4]
    M = str(int(whenstr[4:6]))
    D = str(int(whenstr[6:8]))
    h = str(int(whenstr[8:10]))
    m = str(int(whenstr[10:12]))
    s = str(int(whenstr[12:14]))
    timestamp = Y+'.'+M+'.'+D+'_'+h+'.'+m+'.'+s

    return timestamp


##############################################################################
##############################################################################



def whenstr_from_tstamp(tstamp):
    ymdhms = tstamp.replace('_',' ').replace('.',' ').split()
    if len(ymdhms) != 6:
        raise ValueError('tstamp must be YYYY.M.D_h.m.s -- you entered: '+tstamp)
    whenstr = ymdhms.pop(0)
    if len(whenstr) != 4:
        raise ValueError('Check your timestamp. Year must be 4 digits. You entered: ')
    for x in ymdhms:
        whenstr = whenstr + '%02i' % int(x)
    return whenstr


##############################################################################
##############################################################################


def validate_datestring(datestring):
    """

    Checks whether a datestring is a 4, 6, or 8 character string.
    Checks that the date specified is valid (is a real date) and is in the past.

    """

    if type(datestring) is not str:
        raise TypeError("Your datestring must be a string! You passed in a %s"%type(datestring))

    if len(datestring) > 8:
        raise ValueError("%s is too long. Try YYYYMMDD."%datestring)

    if len(datestring) < 4:
        raise ValueError( "%s is too short. Must be at least YYYY."%datestring)

    if len(datestring) == 5:
        raise ValueError("%s is invalid. You must specify the full month number: YYYYMM"%datestring)

    if len(datestring) == 7:
        raise ValueError("%s is invalid. You must specify the full day number: YYYYMMDD"%datestring)

    if len(datestring) == 8:
        try:
            datetime.strptime(datestring, '%Y%m%d')
        except:
            raise ValueError('%s is not a valid YYYYMMDD date'%datestring)
        if datetime.strptime(datestring, '%Y%m%d') > datetime.utcnow():
            raise ValueError('%s is in the future!'%datestring)

    if len(datestring) == 6:
        try:
            datetime.strptime(datestring, '%Y%m')
        except:
            raise ValueError('%s is not a valid YYYYMM date'%datestring)
        if datetime.strptime(datestring, '%Y%m') > datetime.utcnow():
            raise ValueError('%s is in the future!'%datestring)

    if len(datestring) == 4:
        try:
            datetime.strptime(datestring, '%Y')
        except:
            raise ValueError('%s is not a valid YYYY date'%datestring)
        if datetime.strptime(datestring, '%Y') > datetime.utcnow():
            raise ValueError('%s is in the future!'%datestring)

    return


##############################################################################
##############################################################################


def datestring_to_datedirs(datestring_in, ftype, verbose=True):
    """
    paths = datestring_to_datedirs(datestring, ftype, **verbose)

    Returns a list of strings of the full path to
    Raw or Reduced data directories that match an input date string.
    If verbose=True (default), will ask you to OK if there are multiple directories

    - datestring : yyyy(mmdd) or list of such strings
    - ftype : 'Raw' , 'Reduced'
    - paths : list of path strings

    """

    if 'GPI_TELEM_LOCALDIR' not in os.environ:
            raise Exception('No env variable GPI_TELEM_LOCALDIR. Are you on the right computer?')

    if type(datestring_in) is str:
        datestring_in = [datestring_in]

    if (type(datestring_in) is not list) or (type(datestring_in[0]) is not str):
        raise TypeError('input must be a string or list of strings. You passed a %s of %s'\
        %(type(datestring_in), type(datestring_in[0])))

    datedirs = []

    for datestring in datestring_in:  #
        validate_datestring(datestring)

        if ftype.lower() == 'raw':
            d = os.path.join(os.getenv('GPI_TELEM_DROPBOX'),'Raw')
        elif ftype.lower() == 'reduced':
            d = os.path.join(os.getenv('GPI_TELEM_DROPBOX'),'Reduced')

        thisdatedirs = sorted(glob.glob(os.path.join(d, datestring+'*')))
        print '\n'.join(thisdatedirs)

        if not thisdatedirs:
            raise ValueError('No folders match pattern '+datestring)
        else:
            datedirs += thisdatedirs

    if len(datedirs) > 1 and verbose:
        qstr = 'Your datestring will check '+str(len(datedirs))+' date folders. Do you want to continue? '
        res = query_yes_no(qstr, default='yes')
        if res:
            print 'OK...working through '+str(len(datedirs))+' date folders.\n'
        else:
            sys.exit(0)

    return sorted(datedirs)



##############################################################################
##############################################################################



def check_info_vs_filenames(datestring, verb=True):
    """
    For a given date or date range, finds all "info" files, and checks
    to make sure that all files named in the "info" file header are
    present in the folder. If they're all present, then the set is ready
    to upload to the DB.

    Writes out either "raw_downloaded.txt" or "raw_not_downloaded.txt"
    files to database summary directory.

    Returns a dictionary of missing ugp files for each date; dictionary is
    empty if all data listed in the current info files is present.

    Caveat: obviously, if an "info" file is not yet downloaded, this program
    doesn't know about the entire telem set, and will not report it as missing.

    INPUT:
    - datestring : string or list of strings specifying dates to check.
            Format of strings is YYYYMMDD.
    """


    datedirs = datestring_to_datedirs(datestring,'Raw',verbose=verb)

    to_reduce = {}

    for datedir in sorted(datedirs):
        datestring = datedir.split('/')[-1]
        filelist = sorted(glob.glob(os.path.join(datedir,'*_info.fits')))

        if len(filelist) == 0:
            print 'No aoraw files found'
            return ([],[])

        ugps_missing = {}
        ugps_in_header = {}
        ugps_present = {}

        print "Starting raw file check for "+datestring+"...."

        # for each telem set, make sure all ugp products got sent to AORAW DB
        for f in filelist:

            fshort = f.split('/')[-1] # ugp_When...info.fits
            tstamp = fshort[9:-10] # 201..._...' timestamp

            #print tstamp

            hdulist = fits.open(f)
            h = hdulist[0].header
            ugps_header = h['DTYPE*'] # ugp product files

            #print '------ This is what the header has:'
            #print ugps_header.values() # print

            # make sure all the DTYPE names from the header also show up in the DB query
            badlist = []
            typelist = []
            goodlist = []

            for a_head_ugp in ugps_header.values():
                exists = glob.glob(os.path.join(datedir,a_head_ugp+'.fits'))
                ugpfiletype = a_head_ugp.split(tstamp+'_')[-1]
                typelist.append(ugpfiletype)
                if exists==[]:
                    #print 'bad tstamp ' + tstamp
                    #print 'bad file '+ugpfiletype+'\n'
                    badlist.append( ugpfiletype )
                else:
                    goodlist.append( ugpfiletype )

            if badlist:
                ugps_missing[ tstamp ] = badlist

            if goodlist:
                ugps_present[ tstamp ] = goodlist

            ugps_in_header[ tstamp ] = typelist



        # now save summary files
        summary_dir = os.path.join(os.environ['GPI_TELEM_LOCALDIR'],'database', datestring)
        if not os.path.isdir(summary_dir):
            os.mkdir(summary_dir)


        name_present = 'raw_downloaded.txt'
        name_missing = 'raw_not_downloaded.txt'

        fname0 = os.path.join(summary_dir,name_present)
        if ugps_present.keys():
            print 'Saving summary file: '+fname0+' \n'
            f = open(fname0,'w')
            f.write( "Info DTYPES match file list for these timestamps\n" )
            f.write( '\nChecked on ' + time.ctime() + '\n\n')
            f.write( '\n'.join(ugps_present.keys()) )
            f.close()
        else:
            try: # remove name_present files from previous checks, if present
                os.remove( fname0 )
            except:
                pass


        fname1 = os.path.join(summary_dir,name_missing)
        if ugps_missing.keys():
            print "Some files not yet downloaded from Dropbox"
            # save a summary of missing files
            print 'Saving summary file: '+fname1+'\n'
            f = open(fname1,'w')
            f.write("Some DTYPES from info are not here.\n")
            f.write( '\nChecked on  ' + time.ctime() + '\n\n')
            f.write( '\n'.join(ugps_missing.keys()) )
            f.close()

        else:
            print "Hurray: all DTYPES listed in present info files exist.\n"
            try: # remove missing_raw files from previous checks, if present
                os.remove( fname1 )
            except:
                pass



    return ugps_missing



##############################################################################
##############################################################################



def check_filesizes(date):
    """
    This is a just a starter sketch.
    Don't know what sizes of all the files should be, yet. So don't use.
    """

    #date_dir = os.path.join('/', 'home','sda', 'Dropbox (GPI)','GPI_Telemetry','Raw', date)
    date_dir =  os.path.join(os.getenv('GPI_CODE'),'gpilib','aotelem','Raw',date) # symlink
    filelist = sorted(glob.glob(os.path.join(date_dir,'*_info.fits')))

    if len(filelist) == 0:
        print 'No aoraw files found'
        return ([],[],[])

    timestamps_missing_db = []
    ugps_missing = {}
    ugps_in_header = {}
    ugps_present = {}

    print "Starting file check...."

    # for each telem set, make sure all ugp products got sent to AORAW DB
    for f in filelist:
        #print '\n\n***** Starting file: ' + f + '\n'

        fshort = f.split('/')[-1] # ugp_When...info.fits
        tstamp = fshort[9:-10] # 201..._...' timestamp

        hdulist = fits.open(f)
        h = hdulist[0].header
        ugps_header = h['DTYPE*'] # ugp product files

        #print '------ This is what the header has:'
        #print ugps_header.values() # print

        # make sure all the DTYPE names from the header also show up in the DB query
        badlist = []
        typelist = []
        goodlist = []

        for a_head_ugp in ugps_header.values():
            fname = os.path.join(date_dir, a_head_ugp+'.fits')
            #print fname
            hdulist2 = fits.open(fname)
            h2 = hdulist2[0].header
            s = os.stat(fname)
            fsize = float(s.st_size)
            #print 'fsize = '+str(fsize)
            for ct in np.arange(h2['NAXIS'])+1:
                nkey = 'NAXIS'+str(ct)
                #print '  / '+str(h2[nkey])
                fsize = fsize /  h2[nkey]
                #print 'fsize = '+str(fsize)
            #print '  / '+str( abs(h2['BITPIX']) ) + '  / 8'
            fsize = fsize / ( abs(h2['BITPIX']) / 8.)
            if fsize > 1.2:
                print fname
                print fsize

        return



##############################################################################
##############################################################################



def check_info_vs_db(datestring,db, verb=True):

    """
    For a given date or date range, finds all "info" files, and checks
    to make sure that all files named in the "info" file header are
    present in the DB.

    Writes out either "raw_in_db.txt" or "raw_not_in_db.txt"
    files to database summary directory.

    Returns a dictionary of missing ugp files for each date; dictionary is
    empty if all data listed in the current info files is present.

    Caveat: obviously, if an "info" file is not yet downloaded, this program
    doesn't know about the entire telem set, and will not report it as missing.
    """

    datedirs = datestring_to_datedirs(datestring,'Raw',verbose=verb)

    to_reduce = {}

    for datedir in sorted(datedirs):
        datestring = datedir.split('/')[-1]
        print "\n---------------------------------------------\n"

        # local directory for saving summaries
        summary_dir = os.path.join(os.environ['GPI_TELEM_LOCALDIR'],'database', datestring)
        if verb:
            print 'summary dir is: '+summary_dir
        if not os.path.isdir(summary_dir):
            os.mkdir(summary_dir)

        filelist = sorted(glob.glob(os.path.join(datedir,'*_info.fits')))

        timestamps_missing_db = []
        ugps_missing_db = {}
        ugps_in_header = {}
        ugps_in_db = {}

        if verb:
            print "Starting file check...."

        # for each telem set, make sure all ugp products got sent to AORAW DB
        for f in filelist:
            #print '\n\n***** Starting file: ' + f + '\n'

            fshort = f.split('/')[-1] # ugp_When...info.fits
            tstamp = fshort[9:-10] # 201..._...' timestamp

            hdulist = fits.open(f)
            h = hdulist[0].header
            ugps_header = h['DTYPE*'] # ugp product files

            #print '------ This is what the header has:'
            #print ugps_header.values() # print

            # query all ugp product files in DB for that timestamp.
            prefix = """SELECT DATAFILE from AORAW_FILES
            where DATAFILE like "%"""
            dbquery = prefix + tstamp  + "%\""

            try:
                dbres = db._do_query(dbquery)
                ugps_db = dbres['DATAFILE']
                #print " ----- This is what the DB returns"
                #print ugps_db
            except:
                if verb:
                    print '*********************************'
                    print "No ugp in DB: "+tstamp
                    print '*********************************\n'
                ugps_db = []
                timestamps_missing_db.append(tstamp)

            # make sure all the DTYPE names from the header also show up in the DB query
            badlist = []
            typelist = []
            goodlist = []
            if len(ugps_header.values()) != len(ugps_db):
                for a_head_ugp in ugps_header.values():
                    ugpfiletype = a_head_ugp.split(tstamp+'_')[-1]
                    typelist.append(ugpfiletype)
                    if a_head_ugp not in ugps_db:
                        #print 'bad tstamp ' + tstamp
                        #print 'bad file '+ugpfiletype+'\n'
                        badlist.append( ugpfiletype )
                    else:
                        goodlist.append( ugpfiletype )

            if badlist:
                ugps_missing_db[ tstamp ] = badlist
                ugps_in_header[ tstamp ] = typelist
                ugps_in_db[ tstamp ] = goodlist


        if verb:
            if len(timestamps_missing_db):
                print "Boo - Entire timestamps missing from AO_FILES DB"
                print timestamps_missing_db
                print '\n'
            else:
                print "Hurray: All info file timestamps are in DB!"

        if len(ugps_missing_db):
            if verb:
                print "Boo - UGP files missing from AO_FILES DB:"
            for key in ugps_missing_db.keys():
                to_reduce[key] = sorted(ugps_missing_db[ key ])
                if verb:
                    print key+':  missing keys vs. good keys'
                    print sorted(ugps_missing_db[ key ])
                    print sorted(ugps_in_db[ key ])
                    print '\n'
        else:
            if verb:
                print "Hurray: No ugp files are missing from AO_FILES DB!"


        # write out summary files
        fname_good = os.path.join(summary_dir,'raw_in_db.txt')
        fname_bad = os.path.join(summary_dir,'raw_not_in_db.txt')

        if not timestamps_missing_db and not ugps_missing_db:
            print 'Saving summary file: '+ fname_good +'\n'
            f = open( fname_good,'w' )
            f.write( "All timestamps are present in AO_FILES DB\n\n" )
            f.write( "All individual UGP files are in AO_FILES DB\n" )
            f.write( '\nChecked on ' + time.ctime() + '\n')
            f.close()
            try: # remove raw_not_in_db files from previous checks, if present
                os.remove( fname_bad )
            except:
                pass

        else:
            print 'Saving summary file: '+ fname_bad +'\n'
            f = open( fname_bad ,'w')

            f.write("Entire timestamps missing from AO_FILES DB\n")
            f.write( str(timestamps_missing_db) + '\n\n')

            f.write( "\n\nIndividual UGP files missing from AO_FILES DB\n\n" )
            for key in ugps_missing_db.keys():
                f.write( key +'\n')
                f.write('-- DB missing dtypes --\n' )
                f.write( str(sorted(ugps_missing_db[ key ])) + '\n' )
                f.write('-- DB present dtypes --\n' )
                f.write( str(sorted(ugps_in_db[ key ])) + '\n\n' )

            f.write( '\nChecked on  ' + time.ctime() + '\n')
            f.close()


    return to_reduce



##############################################################################
##############################################################################



def find_raw_to_reduce(datestring, db, verb=True):
    """
    only create error budgets for closed loop on-sky files
    ie: where we're not using the ASU & not @ Zenith

    to_red_dict = find_raw_to_reduce(datestring,db)

    INPUT
    - datestring : yyyymm(dd)
    - db : gpifilesdb instance

    OUTPUT
    - to_red_dict: dictionary of data not in AORED database
        keys = YYYYMMDD
        values = list of [ ugp_..._phase, ARTSRC, OBJNAME, AOCLOOP ] for all missing datasets
    """

    datedirs = datestring_to_datedirs(datestring,'Raw', verbose=verb)

    to_reduce = {}

    for datedir in sorted(datedirs):
        datestring = datedir.split('/')[-1]
        print '----------------------------------------------------------\n'
        print datestring

        prefix = """
        SELECT WHENSTR, ARTSRC, OBJNAME from AORAW
        where WHENSTR like '"""
        query = prefix + datestring + "%'\n order by WHENSTR"

        db_aorawlist = db._do_query(query)
        if db_aorawlist is None:
            raise Exception('There are no AORAW files in the DB for date: '+datestring)


        # make a list of all the files that should be there but *aren't*
        prefix = """
        SELECT AORAW_FILES.DATAFILE, ARTSRC, OBJNAME, AOCLOOP from AORAW
        left join AORAW_FILES on AORAW.UID = AORAW_FILES.RAWID
        where AORAW.UID not in (SELECT RAWID from AORAW2RED) and
        OBJNAME not like 'Zenith' and
        ARTSRC = 0 and
        AOCLOOP = 1 and
        DATAFILE like "%_phase" and
        WHENSTR like '"""
        query = prefix + datestring +"%'\n order by WHENSTR"

        db_toreduce = db._do_query(query)



        # make a list of all the files meeting the criteria to be in the DB
        prefix = """
        SELECT AORAW_FILES.DATAFILE, ARTSRC, OBJNAME, AOCLOOP from AORAW
        left join AORAW_FILES on AORAW.UID = AORAW_FILES.RAWID
        where OBJNAME not like 'Zenith' and
        ARTSRC = 0 and
        AOCLOOP = 1 and
        DATAFILE like "%_phase" and
        WHENSTR like '"""
        query = prefix + datestring +"%'\n order by WHENSTR"

        db_qualified = db._do_query(query)



        summary_dir = os.path.join(os.environ['GPI_TELEM_LOCALDIR'],'database', datestring)
        if not os.path.isdir(summary_dir):
            os.mkdir(summary_dir)


        donename = 'aored_created_and_in_db.txt'
        todoname = 'aored_not_in_db.txt'
        qualifiedname = 'telems_qualified_for_db.txt'

        if db_toreduce is None:
            fname = os.path.join(summary_dir,donename)
            print 'All qualifying AORAW files in DB have AORED DB entries, too.'
            print 'Saving summary: '+fname+' \n'

            f = open(fname,'w')
            f.write( '# All qualifying AORAW files in DB have AORED DB entries too \n' )
            f.write( '# Checked on ' + time.ctime() + '\n' )
            f.close()

            try:
                os.remove( os.path.join(summary_dir,todoname) )
            except:
                pass

        else:
            fname = os.path.join(summary_dir,todoname)
            print 'Saving summary file of ugp that need processed to aored: '+fname+' \n'

            f = open(fname,'w')
            f.write( '# These ao telem sets need to be reduced with IDL errorbudget and uploaded to DB \n' )
            f.write( '# Checked on ' + time.ctime() + '\n' )
            for ugp_phase in db_toreduce['DATAFILE']:
                blah = ugp_phase.split('_')
                tstamp = blah[2]+'_'+blah[3]
                f.write( tstamp+'\n')
            f.close()

            to_reduce[datestring] = db_toreduce

            try:
                os.remove( os.path.join(summary_dir,donename) )
            except:
                pass

        # save a file with all the file names qualified to be in the DB
        try:
            os.remove( os.path.join(summary_dir,qualifiedname) )
        except:
            pass
        fname = os.path.join(summary_dir,qualifiedname)
        f = open(fname,'w')
        if db_qualified is None:
            f.write('# No files qualify for the DB\n')
            f.write( '# Checked on ' + time.ctime() + '\n' )
        else:
            f.write( '# The following files meet the criteria to be included in the DB \n' )
            f.write( '# Checked on ' + time.ctime() + '\n' )
            for ugp_phase in db_qualified['DATAFILE']:
                blah = ugp_phase.split('_')
                tstamp = blah[2]+'_'+blah[3]
                f.write( tstamp+'\n')
        f.close()

    return to_reduce



##############################################################################
##############################################################################



def push_aoraw_2db(datestring, db, verb=True):
    """
    For a given date (range), pushes all the ugp files in that directory to the database.
    Returns a dictionary of files that failed, if applicable

    CAREFUL: This assumes that the raw_downloaded.txt summary file is correct!

    push_aoraw_2db(datestring,db)
    - datestring : yyyymm(dd)
    - db : gpifilesdb instance
   """

    datedirs = datestring_to_datedirs(datestring,'Raw', verbose=verb)

    allbadfiles = {}
    for datedir in datedirs:
        thisdatestring = datedir.split('/')[-1]

        # first sanity check to make sure you've checked whether files are all there
        summary_dir = os.path.join(os.environ['GPI_TELEM_LOCALDIR'],'database',thisdatestring)
        if ( not os.path.isfile(os.path.join(summary_dir,'raw_downloaded.txt')) ) or\
           os.path.isfile(os.path.join(summary_dir, 'raw_not_downloaded.txt')):
            print 'WARNING: skipping '+thisdatestring+' b/c no raw_downloaded summary file'+\
              ' or b/c raw_not_downloaded summary file exists.'
            allbadfiles[thisdatestring] = "Didn't try b/c some raw data missing."
            sys.exit(-1)


        # "info" files are what's use to upload raw to DB
        filelist = sorted(glob.glob(os.path.join(datedir,'*_info.fits')))

        if len(filelist) == 0:
            print 'No aoraw files found for '+thisdatestring
            continue

        badfiles = []
        for f in filelist:
            fshort = f.split('/')[-1]
            try:
                db.updateRawAOProduct(f)
                print "Finished "+fshort
            except:
                print '*********************************'
                print "Error uploading: "+f
                print '*********************************\n'
                badfiles.append(f)

        if badfiles: allbadfiles[thisdatestring] = badfiles

    if allbadfiles:
        print 'some files could not be uploaded'
        print allbadfiles
    else:
        print 'Hurray! No AORAW uploads threw an error.'

    return allbadfiles



##############################################################################
##############################################################################


def push_aoraw_2db_update(tstamps, db, verb=True):
    """
    badtimestamps = push_aoraw_2db_update(tstamps, db)

    Sometimes only some files from AO raw data upload properly.
    This takes a (list of) YYYY.M.D_h.m.s and
    tries to re-upload any missing raw (ugp) data to the DB.
    Outputs a list of timestamps that still didn't work.

    - tstamps : (list of) YYYY.M.D_h.m.s to try to fix
    - db : GPIFilesDB instance
    - badtimestamps : list of YYYY.M.D_h.m.s that refuse to be fixed

    """

    # if single ts input, make it a list for the loop
    if isinstance(tstamps,str): tstamps = [tstamps]

    badts = []

    for tstamp in tstamps:
        ws = whenstr_from_tstamp(tstamp)
        datestr = ws[:8]
        datedirs = datestring_to_datedirs(datestr, 'Raw', verbose=verb)
        if len(datedirs) !=1:
            sys.exit("You shouldn't be running multiple directories here!")
        datedir = datedirs[0]

        filelist = sorted(glob.glob(os.path.join(datedir,'ugp_When_'+tstamp+'_info.fits')))
        if len(filelist) == 0:
            sys.exit("no files found for "+tstamp)
        if len(filelist) > 1:
            sys.exit("multiple 'info' files with the same name?!?!")

        try:
            db.updateRawAOProduct_addmissingdtypes(filelist[0])
            print '*********************************'
            print "Re-uploaded  " + tstamp
            print '*********************************\n'
        except:
            print '*********************************'
            print "Error uploading: " + tstamp
            print '*********************************\n'
            badts.append(tstamp)

    return badts



##############################################################################
##############################################################################



def remove_old_aored_from_DB(datestring, db, verb=True):
    """
    For a given date (range), remove any AORED and AORAW2RED entries
    that don't match the UID of what is currently saved in Dropbox.
    ie: remove entries from old data reductions

    remove_old_aored_from_DB(datestring, db, [verb=True/False])
    - datestring : yyyy(mmdd)
    - db : gpifiles db instance
    - verb : verbose output? default True
    """
    import hashlib

    datedirs = datestring_to_datedirs(datestring,'Reduced', verbose=verb)

    allbadfiles = {}
    for datedir in datedirs:
        thisdatestring = datedir.split('/')[-1]

        d = os.path.join(os.path.join(os.getenv('GPI_TELEM_DROPBOX'),'Reduced', thisdatestring))
        filelist = sorted(glob.glob(os.path.join(d,'*_lap_eb.fits')))
        if len(filelist) == 0:
            print 'No aored files found for '+thisdatestring
            allbadfiles[thisdatestring] = 'No AORED files at all'
            continue

        badfiles = []
        for f in filelist:
            fshort = f.split('/')[-1]
            tstmp = fshort[11:-12]

            hdulist = fits.open(f)
            uid = hashlib.md5(open(hdulist.filename(), 'rb').read()).hexdigest()

            q1 = "SELECT UID from AORED where datafile like '%" + tstmp + "%'"
            q2 = " and UID not like '" +  str(uid) + "'"
            uids_other = db._do_query(q1+q2)

            if uids_other is None:
                print 'no duplicate for '+tstmp
                continue

            for olduid in uids_other[0:1].UID:
                try:
                    db._do_delete("""DELETE FROM AORAW2RED Where redid='%s'""" % (olduid))
                    db._do_delete("""DELETE FROM AORED WHERE uid='%s'""" % (olduid))
                    print "old entry removed for '%s'" % (tstmp)
                except:
                    print '\n*********************************'
                    print "Error removing duplicate for: "+tstmp
                    print '*********************************\n'
                    badfiles.append((tstmp,olduid))

        if badfiles: allbadfiles[thisdatestring] = badfiles

    if allbadfiles:
        print 'some AORED files threw errors when removing:'
        print allbadfiles
    else:
        print 'Hurray! Old AORED removed successfully.'

    return allbadfiles



##############################################################################
##############################################################################



def push_aored_2db(datestring, db, verb=True):
    """
    For a given date (range), pushes all the DB-qualified AORED files in that directory to the DB.
    Returns a dictionary of files that failed, if applicable

    push_aored_2db(datestring,db)
    - datestring : yyyymm(dd)
    - db : gpifilesdb instance
   """

    datedirs = datestring_to_datedirs(datestring,'Reduced', verbose=verb)

    allbadfiles = {}
    for datedir in datedirs:
        thisdatestring = datedir.split('/')[-1]

        d = os.path.join(os.path.join(os.getenv('GPI_TELEM_DROPBOX'),'Reduced', thisdatestring))
        filelist = sorted(glob.glob(os.path.join(d,'*_lap_eb.fits')))
        if len(filelist) == 0:
            print 'No aored files found for '+thisdatestring
            allbadfiles[thisdatestring] = 'No AORED files at all'
            continue


        # make a list of qualified files
        summary_dir = os.path.join(os.environ['GPI_TELEM_LOCALDIR'],'database',thisdatestring)
        qual_file = os.path.join(summary_dir, 'telems_qualified_for_db.txt')
        if not os.path.isfile(qual_file):
            print 'WARNING: no telems_qualified_for_db.txt in '+thisdatestring+'... skipping'
            continue
        f = open(qual_file,'r')
        tmp = f.read()
        qual_list = tmp.split('\n')
        qual_list = [x for x in qual_list if not x.startswith('#') ]
        qual_list = filter(None,qual_list) #remove blank entries

        badfiles = []
        for f in filelist:
            fshort = f.split('/')[-1]
            tstmp = fshort[11:-12]
            if not tstmp in qual_list:
                print '\n*********************************\n'
                print '------- '+tstmp+' not qualified for DB. Skipping...'
                print '*********************************\n'
                continue

            try:
                # use replace=1 if you want to overwrite existing DB entries for these files
                db.updateReducedAOProduct(f)
                #print "Not actually uploading aored to db!"
                print tstmp+': OK '
            except:
                print '\n*********************************'
                print "Error uploading: "+f
                print '*********************************\n'
                badfiles.append(f)

        if badfiles: allbadfiles[thisdatestring] = badfiles

        # check that they're uploaded, and update summary file
        res=find_raw_to_reduce(thisdatestring, db)


    if allbadfiles:
        print 'some qualified AORED files were not uploaded'
        print allbadfiles
    else:
        print 'Hurray! No AORED uploads threw an error.'

    return allbadfiles



##############################################################################
##############################################################################


def main():
    # nothing to do here, since this is just a toolkit
    print "I don't do anything by myself"


if __name__ == '__main__':
    main()
