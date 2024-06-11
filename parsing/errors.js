


class ErrorTracker {

    constructor(errorDivID, successMsg, callback=null, successCallback=null) {
        this.hasError = false;
        this.message = null;
        this.successMsg = successMsg;
        this.target = errorDivID;
        this.callback = callback;
        this.successCallback = successCallback;
    }

    setTarget(errorDivID) {
        this.target = errorDivID;
    }

    setCallback(callback) {
        this.callback = callback;
    }

    setSuccessCallback(successCallback) {
        this.successCallback = successCallback;
    }

    defaultCallback(target) {
        if (target) {
            target.innerText = this.message;
            target.style.color = "red";
            target.style.display = "block";
        }
    }

    defaultSuccessCallback(target) {
        if (target) {
            target.innerText = this.message;
            target.style.color = "green";
            target.style.display = "block";
        }
    }

    error(message) {
        this.hasError = true;
        this.message = message;

        const target = document.querySelector(`#${this.target}`);

        if (this.callback === null) {
            this.defaultCallback(target);
        } else {
            this.callback(this.message, target);
        }
    }

    clear() {
        this.hasError = false;
        this.message = this.successMsg;

        const target = document.querySelector(`#${this.target}`);
        
        if (this.successCallback === null) {
            this.defaultCallback(target);
        } else {
            this.successCallback(this.message, target);
        }
    }

}


const tracker = new ErrorTracker(null, "Parsing successful!");