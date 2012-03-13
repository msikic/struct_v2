# include "struct.h"
char single_letter ( char code[]){

    switch ( code[0] ) {
    case 'A':
	switch ( code [1]) {
	case 'L':
	    return 'A';
	    break;
	case 'R':
	    return 'R';
	    break;
	case 'S':
	    switch ( code[2] ) {
	    case 'N':
		return 'N';
		break;
	    case 'P':
		return  'D';
		break;
	    }
	    break;
	}
	break;
    case 'C':
	return 'C'; 
	break;
    case 'G':
	/* the second letter is always L */ 
	switch ( code[2] ) {
	case 'U':
	    return 'E';
	    break;
	case 'N':
	    return  'Q';
	    break;
	case 'Y':
	    return 'G';
	    break;
	}
	break;
    case 'H':
	return  'H';
	break;
    case 'I':
	return  'I';
	break;
    case 'L':
	switch ( code [1]) {
	case 'E':
	    return 'L';
	    break;
	case 'Y':
	    return 'K';
	    break;
	}
	break;
    case 'M':
	return 'M';
	break;
    case 'P':
	switch ( code [1]) {
	case 'H':
	    return 'F';
	    break;
	case 'R':
	    return 'P';
	    break;
	}
	break;
    case 'S':
	return 'S';
	break;
    case 'T':
	switch ( code [1]) {
	case 'H':
	    return 'T';
	    break;
	case 'R':
	    return 'W';
	    break;
	case 'Y':
	    return 'Y';
	    break;
	}
	break;
    case 'V':
	return 'V';
	break;
	
    }


    //fprintf (stdout, "Unrecognized amino acid code: %s.\n", code);
    return 0;
}

