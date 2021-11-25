{
	//int run;
	//std::cout << "\nEnter run number: ";
	//std::cin >> run;



	TRIFIC->cd();
	trifEdE->GetXaxis()->SetTickLength(2);
	trifEdE->GetYaxis()->SetTickLength(2);
	trifEdE->GetXaxis()->SetRangeUser(2000,7000);
	trifEdE->GetYaxis()->SetRangeUser(1000,15000);
	//trifEdE->SetTitle(Form("PID Run %i",(int)run));
	trifEdE->Draw("colz");

}
