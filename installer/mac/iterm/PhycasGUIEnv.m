// -*- mode:objc -*-
// $Id: PhycasGUIApplication.m,v 1.10 2006/11/07 08:03:08 yfabian Exp $
//
/*
 **  PhycasGUIApplication.m
 **
 **  Copyright (c) 2002-2004
 **
 **  Author: Ujwal S. Setlur
 **
 **  Project: iTerm
 **
 **  Description: overrides sendEvent: so that key mappings with command mask  
 **				  are handled properly.
 **
 **  This program is free software; you can redistribute it and/or modify
 **  it under the terms of the GNU General Public License as published by
 **  the Free Software Foundation; either version 2 of the License, or
 **  (at your option) any later version.
 **
 **  This program is distributed in the hope that it will be useful,
 **  but WITHOUT ANY WARRANTY; without even the implied warranty of
 **  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **  GNU General Public License for more details.
 **
 **  You should have received a copy of the GNU General Public License
 **  along with this program; if not, write to the Free Software
 **  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#import "PhycasGUIEnv.h"

@implementation PhycasGUIEnv
 + (id) sharedInstance 
{
	static PhycasGUIEnv * shared = nil;
	if ( !shared ) 
	{
		shared = [[self alloc] init];
	}
	return shared; 
}

- (NSString *) guiArgvZero
{
	return phycasGUIArgvZero;
}
- (void)setGuiArgvZero:(NSString *)newArgvZ
{
	[phycasGUIArgvZero release];
	phycasGUIArgvZero = newArgvZ;
}

- init 
{ 
	return[self initWithName:"default"]; 
} 

- initWithName:(char *)string 
{ 
	self = [super init];
	if ( self ) {
		phycasGUIArgvZero = [NSString stringWithFormat: @"%s", string];
	} 
	return self; 
}
@end
