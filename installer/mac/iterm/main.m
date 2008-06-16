// -*- mode:objc -*-
// $Id: main.m,v 1.1.1.1 2002/11/26 04:56:45 ujwal Exp $
//
//  main.m
//  JTerminal
//
//  Created by kuma on Thu Nov 22 2001.
//  Copyright (c) 2001 Kiichi Kusama. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import "PhycasGUIEnv.h"

int main(int argc, const char *argv[])
{
	NSString *s = [NSString stringWithFormat: @"%s", argv[0]];
	[[PhycasGUIEnv sharedInstance] setGuiArgvZero: s];
    return NSApplicationMain(argc, argv);
}
