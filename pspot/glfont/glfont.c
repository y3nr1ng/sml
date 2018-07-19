//*********************************************************
//GLFONT.CPP -- glFont routines
//Copyright (c) 1998 Brad Fish
//See glFont.txt for terms of use
//November 10, 1998
//*********************************************************

// #include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif
#include "glfont.h"

//*********************************************************
//Variables
//*********************************************************

//Current font
GLFONT *glFont;

//*********************************************************
//Functions
//*********************************************************
int glFontCreate (GLFONT *Font, char *FileName, int Tex)
{
	FILE          *Input;
	unsigned char *TexBytes;
	int Num;

	//Open font file
	if ((Input = fopen(FileName, "rb")) == NULL) {
		printf("!!! glFontCreate: cannot open file: %s\n", FileName);
		return FALSE;
	}

	//Read glFont structure
	fread(Font, 20, 1, Input);
	fseek(Input, 4, SEEK_CUR);

	//Save texture number
	Font->Tex = Tex;

	//Allocate memory for characters
	Num = Font->IntEnd - Font->IntStart + 1;
	if ((Font->Char = (GLFONTCHAR *)malloc(
		sizeof(GLFONTCHAR) * Num)) == NULL)
		return FALSE;

	//Read glFont characters
	fread(Font->Char, sizeof(GLFONTCHAR), Num, Input);

	//Allocate memory for texture data for the font file
	Num = Font->TexWidth * Font->TexHeight * 2;
	if ((TexBytes = malloc(Num)) == NULL)
		return FALSE;
	
	//Read texture data
	fread(TexBytes, sizeof(char), Num, Input);

	//Set texture attributes
	glBindTexture(GL_TEXTURE_2D, Font->Tex);  
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);  
	
	//Create texture
	glTexImage2D(GL_TEXTURE_2D, 0, 2, Font->TexWidth,
		Font->TexHeight, 0, GL_LUMINANCE_ALPHA, 
		GL_UNSIGNED_BYTE, (void *)TexBytes);

	//Clean up
	free(TexBytes);
	fclose(Input);

	//Return pointer to new font
	return TRUE;
}
//*********************************************************
void glFontDestroy (GLFONT *Font)
{
	//Free character memory
	free(Font->Char);
}
//*********************************************************
void glFontBegin (GLFONT *Font)
{
	//Save pointer to font structure
	if (Font->Char != NULL)
		glFont = Font;
	else
		glFont = NULL;
	
	//Bind to font texture
	glBindTexture(GL_TEXTURE_2D, Font->Tex);
}
//*********************************************************
void glFontEnd (void)
{
	//Font no longer current
	glFont = NULL;
}
//*********************************************************
void glFontTextOut (char *String, float x, float y, float z)
{
	int Length, i;
	GLFONTCHAR *Char;

	//Return if we don't have a valid glFont 
	if (glFont == NULL)
		return;
	
	//Get length of string
	Length = strlen(String);
	
	//Begin rendering quads
	glBegin(GL_QUADS);
	
	//Loop through characters
	for (i = 0; i < Length; i++)
	{
		//Get pointer to glFont character
		Char = &glFont->Char[(int)String[i] -
			glFont->IntStart];
		
		//Specify vertices and texture coordinates
		glTexCoord2f(Char->tx1, Char->ty1);
		glVertex3f(x, y, z);
		glTexCoord2f(Char->tx1, Char->ty2);
		glVertex3f(x, y - Char->dy, z);
		glTexCoord2f(Char->tx2, Char->ty2);
		glVertex3f(x + Char->dx, y - Char->dy, z);
		glTexCoord2f(Char->tx2, Char->ty1);
		glVertex3f(x + Char->dx, y, z);
	
		//Move to next character
		x += Char->dx;
	}

	//Stop rendering quads
	glEnd();
}
//*********************************************************

//End of file



