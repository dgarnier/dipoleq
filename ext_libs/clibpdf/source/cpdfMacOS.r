/* Apparently, Apple event apps must have a SIZE resource */
/* Do we really need this for MacOS 8/7?  What's a Rez ? */

resource 'SIZE' (-1, purgeable) {
                reserved,
                acceptSuspendResumeEvents,
                reserved,
                canBackground,
                multiFinderAware,
                backgroundOnly,
                dontGetFrontClicks,
                ignoreChildDiedEvents,
                is32BitCompatible,
                isHighLevelEventAware,
                localAndRemoteHLEvents,
                isStationeryAware,
                dontUseTextEditServices,
                reserved,
                reserved,
                reserved,
                524288,
                524288
};

	/* backgroundOnly was backgroundAndForeground */

