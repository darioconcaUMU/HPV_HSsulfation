function [tracks, track_props, discarded, discarded_props] = filter_duration(tracks, track_props, min_duration)

short_tracks = [track_props(:).duration] < min_duration;

discarded = tracks(short_tracks);
discarded_props = track_props(short_tracks);

tracks(short_tracks) = [];
track_props(short_tracks) = [];
